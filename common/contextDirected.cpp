#include "contextDirected.h"
#include <boost/graph/graphml.hpp>
#include <fstream>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/scoped_array.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/serialization/map.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
namespace multistateTurnip
{
	namespace ContextDirectedImpl
	{
		template<typename Key, ContextDirected::internalGraph::vertices_size_type Ret> class constant_property_map_vertices_size_type : public boost::put_get_helper<ContextDirected::internalGraph::vertices_size_type, constant_property_map_vertices_size_type<Key, Ret> > 
		{
		public:
			typedef Key key_type;
			typedef ContextDirected::internalGraph::vertices_size_type reference;
			typedef ContextDirected::internalGraph::vertices_size_type value_type;

			typedef boost::readable_property_map_tag category;

			constant_property_map_vertices_size_type(){}

			reference operator[](const Key&) const 
			{
				return Ret;
			}
		};
		template<typename Key, int Ret> class constant_property_map_int : public boost::put_get_helper<int, constant_property_map_int<Key, Ret> > 
		{
		public:
			typedef Key key_type;
			typedef int reference;
			typedef int value_type;

			typedef boost::readable_property_map_tag category;

			constant_property_map_int(){}

			reference operator[](const Key&) const 
			{
				return Ret;
			}
		};

		struct twoDArray
		{
			int* base;
			std::size_t dim;
			struct twoDArrayInternal
			{
				twoDArrayInternal(int* base)
					:base(base)
				{};
				int& operator[](std::size_t j)
				{
					return *(base + j);
				}
				const int& operator[](std::size_t j) const
				{
					return *(base + j);
				}
				int* base;
			};
			twoDArrayInternal operator[](std::size_t i) const
			{
				return twoDArrayInternal(base + dim*i);
			};
		};
	}
	ContextDirected& ContextDirected::operator=(ContextDirected&& other)
	{
		graph.swap(other.graph);
		source = other.source;
		sink = other.sink;
		distributions = std::move(other.distributions);

		edgeResidualCapacityVector.swap(other.edgeResidualCapacityVector);
		capacityVector.swap(other.capacityVector);
		vertexPredecessorVector.swap(other.vertexPredecessorVector);
		colorVector.swap(other.colorVector);
		distanceVector.swap(other.distanceVector);

		std::swap(threshold, other.threshold);
		return *this;
	}
	ContextDirected::ContextDirected(ContextDirected&& other)
	{
		graph.swap(other.graph);
		source = other.source;
		sink = other.sink;
		distributions = std::move(other.distributions);

		edgeResidualCapacityVector.swap(other.edgeResidualCapacityVector);
		capacityVector.swap(other.capacityVector);
		vertexPredecessorVector.swap(other.vertexPredecessorVector);
		colorVector.swap(other.colorVector);
		distanceVector.swap(other.distanceVector);

		std::swap(threshold, other.threshold);
	}
	ContextDirected::ContextDirected(boost::shared_ptr<const inputGraph> inputDirectedGraph, int source, int sink, const std::vector<capacityDistribution>& inputDistributions, const mpfr_class& threshold)
		: source(source), sink(sink)
	{
		std::size_t nVertices = boost::num_vertices(*inputDirectedGraph);
		int minVertexIndex = std::min(source, sink);
		int maxVertexIndex = std::max(source, sink);
		if(minVertexIndex < 0 || maxVertexIndex >= (int)nVertices)
		{
			throw std::runtime_error("Input interestVertices was out of range");
		}
		{
			std::size_t nEdges = boost::num_edges(*inputDirectedGraph);
			if(this->distributions.size() != nEdges)
			{
				throw std::runtime_error("Input distributions must have an entry for every edge");
			}
			
			capacityVector.resize(nEdges);
			if(boost::connected_components(*inputDirectedGraph, &(capacityVector[0])) > 1)
			{
				throw std::runtime_error("Input graph must contain a single connected component");
			}
			//Check that the edge indices of the input graph are unique, and are in [0, nEdges)
			std::vector<int> counters(nEdges, 0);
			inputGraph::edge_iterator current, end;
			boost::tie(current, end) = boost::edges(*inputDirectedGraph);
			for(; current != end; current++)
			{
				int edgeIndex = boost::get(boost::edge_index, *inputDirectedGraph, *current);
				if(edgeIndex < 0 || edgeIndex >= (int)nEdges)
				{
					throw std::runtime_error("Edge index was out of range");
				}
				if(counters[edgeIndex] > 0)
				{
					throw std::runtime_error("Edge indices were not unique");
				}
				counters[edgeIndex]++;
			}
		}

		//Construct a new graph, where the edge indices are assigned consecutively
		boost::shared_ptr<internalGraph> internalDirectedGraph(new internalGraph(nVertices));
		inputGraph::edge_iterator start, end;
		boost::tie(start, end) = boost::edges(*inputDirectedGraph);
		int index = 0;

		std::vector<std::pair<mpfr_class, mpfr_class>> alwaysZeroCapacity;
		alwaysZeroCapacity.push_back(std::make_pair(mpfr_class(0), mpfr_class(1)));
		capacityDistribution zeroDistribution(alwaysZeroCapacity);

		for(; start != end; start++)
		{
			//This has all the reverse edges too. 
			int source = boost::source(*start, *inputDirectedGraph), target = boost::target(*start, *inputDirectedGraph);
			std::pair<internalGraph::edge_descriptor, bool> tmp = boost::edge(source, target, *internalDirectedGraph);
			//Doing it this way because I'm not sure if add_edge overwrites the edge_index in the case that the edge already exists. 
			internalGraph::edge_descriptor firstEdge;
			int inputEdgeIndex = boost::get(boost::edge_index, *inputDirectedGraph, *start);
			if(!tmp.second)
			{
				firstEdge = boost::add_edge(source, target, index, *internalDirectedGraph).first;
				//Put in the distribution for the original directed edge
				distributions.push_back(inputDistributions[inputEdgeIndex].makeCopy());
				index++;
			}
			else 
			{
				firstEdge = tmp.first;
				distributions[boost::get(boost::edge_index, *internalDirectedGraph, firstEdge)] = inputDistributions[inputEdgeIndex].makeCopy();
			}

			internalGraph::edge_descriptor secondEdge;
			tmp = boost::edge(target, source, *internalDirectedGraph);
			if(!tmp.second)
			{
				secondEdge = boost::add_edge(target, source, index, *internalDirectedGraph).first;
				distributions.push_back(alwaysZeroCapacity);
				index++;
			}
			else secondEdge = tmp.first;
			//Put in reverse edges
			boost::put(boost::edge_reverse, *internalDirectedGraph, firstEdge, secondEdge);
			boost::put(boost::edge_reverse, *internalDirectedGraph, secondEdge, firstEdge);
		}
		std::size_t nEdgesWithReverse = boost::num_edges(*internalDirectedGraph);
		capacityVector.resize(nEdgesWithReverse);
		edgeResidualCapacityVector.resize(nEdgesWithReverse);
		vertexPredecessorVector.resize(nEdgesWithReverse);
		vertexPredecessorVector.resize(nEdgesWithReverse);
		colorVector.resize(nEdgesWithReverse);
		distanceVector.resize(nEdgesWithReverse);

		this->graph = internalDirectedGraph;
	}
	ContextDirected::ContextDirected()
		:graph(NULL)
	{}
	ContextDirected ContextDirected::emptyContextDirected()
	{
		ContextDirected result;

		boost::shared_ptr<ContextDirected::internalGraph> graph(new ContextDirected::internalGraph(2));
		boost::add_edge(0, 1, 0, *graph);
		boost::add_edge(1, 0, 1, *graph);
		boost::shared_ptr<const ContextDirected::internalGraph> constGraph = boost::static_pointer_cast<const ContextDirected::internalGraph>(graph);
		result.graph.swap(constGraph);

		result.source = 0;
		result.sink = 1;
		
		result.edgeResidualCapacityVector.resize(2);
		result.capacityVector.resize(2);
		result.vertexPredecessorVector.resize(2);
		result.colorVector.resize(2);
		result.distanceVector.resize(2);

		return result;
	}
	ContextDirected ContextDirected::gridContextDirected(int gridDimension, int source, int sink, const std::vector<capacityDistribution>& distributions, const mpfr_class& threshold)
	{
		boost::shared_ptr<ContextDirected::inputGraph> graph(new ContextDirected::inputGraph(gridDimension * gridDimension));
		int counter = 0;
		for(int i = 0; i < gridDimension; i++)
		{
			for(int j = 0; j < gridDimension; j++)
			{
				if(i != gridDimension - 1) 
				{
					boost::add_edge(i + j*gridDimension, i + 1 + j*gridDimension, counter++, *graph);
					boost::add_edge(i + 1 + j*gridDimension, i + j*gridDimension, counter++, *graph);
				}
				if(j != gridDimension - 1)
				{
					boost::add_edge(i + j*gridDimension, i + (j+1)*gridDimension, counter++, *graph);
					boost::add_edge(i + (j+1)*gridDimension, i + j*gridDimension, counter++, *graph);
				}
			}
		}
		return ContextDirected(graph, source, sink, distributions, threshold);
	}
	ContextDirected ContextDirected::completeContextDirected(int nVertices, const std::vector<capacityDistribution>& distributions, const mpfr_class& threshold)
	{
		boost::shared_ptr<ContextDirected::inputGraph> graph(new ContextDirected::inputGraph(nVertices));
		int counter = 0;
		for(int i = 0; i < nVertices; i++)
		{
			for(int j = 0; j < nVertices; j++)
			{
				if(i != j) 
				{
					boost::add_edge(i, j, counter++, *graph);
					boost::add_edge(j, i, counter++, *graph);
				}
			}
		}

		return ContextDirected(graph, 0, 1, distributions, threshold);
	}
	const ContextDirected::internalGraph& ContextDirected::getGraph() const
	{
		return *graph;
	}
	int ContextDirected::getSource() const
	{
		return source;
	}
	int ContextDirected::getSink() const
	{
		return sink;
	}
	ContextDirected ContextDirected::fromFile(std::string path, bool& successful, int source, int sink, std::string& message, const std::vector<capacityDistribution>& distributions, const mpfr_class& threshold)
	{
		std::ifstream input(path);
		if(!input.is_open())
		{
			successful = false;
			return ContextDirected();
		}
		boost::dynamic_properties properties;
		boost::shared_ptr<inputGraph> graph(new inputGraph());

		boost::vector_property_map<int> orderingProperty;
		properties.property("order", orderingProperty);

		boost::vector_property_map<float> xProperty, yProperty;
		properties.property("x", xProperty);
		properties.property("y", yProperty);

		boost::read_graphml(input, *graph, properties);

		successful = true;

		if(std::min(source, sink) < 0 || std::max(source, sink) >= (int)boost::num_vertices(*graph))
		{
			successful = false;
			message = "Invalid vertex indices entered for input interestVertices";
			return ContextDirected();
		}
		
		inputGraph::edge_iterator current, end;
		boost::tie(current, end) = boost::edges(*graph);
		int counter = 0;
		for(; current != end; current++)
		{
			boost::put(boost::edge_index, *graph, *current, counter++);
		}
		return ContextDirected(graph, source, sink, distributions, threshold);
	}
	const capacityDistribution& ContextDirected::getDistribution(int edgeIndex) const
	{
		return distributions[edgeIndex];
	}
	ContextDirected::~ContextDirected()
	{
	}
	std::vector<double>& ContextDirected::getCapacityVector() const
	{
		return capacityVector;
	}
/*	double ContextDirected::getMaxFlow(std::vector<double>& capacities, ContextDirected::internalGraph::vertex_descriptor source, ContextDirected::internalGraph::vertex_descriptor sink) const
	{
		typedef boost::property_map<ContextDirected::internalDirectedGraph, boost::edge_index_t>::const_type edgeIndexMapType;
		typedef boost::property_map<ContextDirected::internalDirectedGraph, boost::vertex_index_t>::const_type vertexIndexMapType;
		typedef boost::iterator_property_map<typename std::vector<double>::iterator, edgeIndexMapType> edgeCapacityMapType;
		typedef boost::iterator_property_map<typename std::vector<ContextDirected::internalDirectedGraph::edge_descriptor>::iterator, vertexIndexMapType> vertexPredecessorMapType;
		typedef boost::iterator_property_map<typename std::vector<boost::default_color_type>::iterator, vertexIndexMapType> colorMapType;
		typedef boost::iterator_property_map<typename std::vector<int>::iterator, vertexIndexMapType> distanceMapType;

		edgeIndexMapType edgeIndexMap = boost::get(boost::edge_index, *directedGraph);
		vertexIndexMapType vertexIndexMap = boost::get(boost::vertex_index, *directedGraph);
		edgeCapacityMapType residualCapacityMap(edgeResidualCapacityVector.begin(), edgeIndexMap);
		edgeCapacityMapType edgeCapacityMap(capacities.begin(), edgeIndexMap);
		vertexPredecessorMapType vertexPredecessorMap(vertexPredecessorVector.begin(), vertexIndexMap);
		colorMapType colorMap(colorVector.begin(), vertexIndexMap);
		distanceMapType distanceMap(distanceVector.begin(), vertexIndexMap);

		return boost::boykov_kolmogorov_max_flow(*directedGraph, source, sink, boost::residual_capacity_map(residualCapacityMap).capacity_map(edgeCapacityMap).predecessor_map(vertexPredecessorMap).color_map(colorMap).distance_map(distanceMap));
	}
	double ContextDirected::getMaxFlow(std::vector<double>& capacities) const
	{
		return getMaxFlow(capacities, (*interestVertices)[0], (*interestVertices)[1]);
	}*/
}
