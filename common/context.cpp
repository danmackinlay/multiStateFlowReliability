#include "context.h"
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
	namespace ContextImpl
	{
		template<typename Key, Context::internalGraph::vertices_size_type Ret> class constant_property_map_vertices_size_type : public boost::put_get_helper<Context::internalGraph::vertices_size_type, constant_property_map_vertices_size_type<Key, Ret> > 
		{
		public:
			typedef Key key_type;
			typedef Context::internalGraph::vertices_size_type reference;
			typedef Context::internalGraph::vertices_size_type value_type;

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
	Context& Context::operator=(Context&& other)
	{
		graph.swap(other.graph);
		interestVertices.swap(other.interestVertices);
		vertexPositions.swap(other.vertexPositions);
		edgeDistances.swap(other.edgeDistances);
		distribution = std::move(other.distribution);

		edgeResidualCapacityVector.swap(other.edgeResidualCapacityVector);
		capacityVector.swap(other.capacityVector);
		directedGraph.swap(other.directedGraph);
		vertexPredecessorVector.swap(other.vertexPredecessorVector);
		colorVector.swap(other.colorVector);
		distanceVector.swap(other.distanceVector);

		std::swap(threshold, other.threshold);
		std::swap(nEdges, other.nEdges);
		return *this;
	}
	Context::Context(Context&& other)
	{
		graph.swap(other.graph);
		interestVertices.swap(other.interestVertices);
		vertexPositions.swap(other.vertexPositions);
		edgeDistances.swap(other.edgeDistances);
		distribution = std::move(other.distribution);

		edgeResidualCapacityVector.swap(other.edgeResidualCapacityVector);
		capacityVector.swap(other.capacityVector);
		directedGraph.swap(other.directedGraph);
		vertexPredecessorVector.swap(other.vertexPredecessorVector);
		colorVector.swap(other.colorVector);
		distanceVector.swap(other.distanceVector);

		std::swap(threshold, other.threshold);
		std::swap(nEdges, other.nEdges);
	}
	Context::Context(boost::shared_ptr<const inputGraph> unorderedGraph, boost::shared_ptr<const std::vector<int> > interestVertices, boost::shared_ptr<std::vector<vertexPosition> > vertexPositions, capacityDistribution&& distribution, const mpfr_class& threshold)
		:interestVertices(interestVertices), vertexPositions(vertexPositions), distribution(std::move(distribution))
	{
		if(interestVertices->size() != 2)
		{
			throw std::runtime_error("Two interest vertices must be entered");
		}

		std::size_t nVertices = boost::num_vertices(*unorderedGraph);
		int minVertexIndex = *std::min_element(interestVertices->begin(), interestVertices->end());
		int maxVertexIndex = *std::max_element(interestVertices->begin(), interestVertices->end());
		if(minVertexIndex < 0 || maxVertexIndex >= nVertices)
		{
			throw std::runtime_error("Input interestVertices was out of range");
		}
		nEdges = boost::num_edges(*unorderedGraph);
		
		edgeResidualCapacityVector.resize(2*nEdges);
		vertexPredecessorVector.resize(2*nEdges);
		capacityVector.resize(2*nEdges, 1);
		vertexPredecessorVector.resize(2*nEdges);
		colorVector.resize(2*nEdges);
		distanceVector.resize(2*nEdges);

		if(boost::connected_components(*unorderedGraph, &(capacityVector[0])) > 1)
		{
			throw std::runtime_error("Input graph must contain a single connected component");
		}
		if(nVertices != vertexPositions->size())
		{
			throw std::runtime_error("Vertex position data had the wrong size");
		}

		boost::shared_ptr<internalGraph> orderedGraph(new internalGraph(nVertices));
		inputGraph::edge_iterator start, end;
		boost::tie(start, end) = boost::edges(*unorderedGraph);
		int index = 0;
		for(; start != end; start++)
		{
			boost::add_edge(start->m_source, start->m_target, index, *orderedGraph);
			index++;
		}

		graph = orderedGraph;
		constructDirectedGraph();
		constructEdgeDistances();
	}
	void Context::constructEdgeDistances()
	{
		const std::size_t nEdges = boost::num_edges(*graph);
		const std::size_t nVertices = boost::num_vertices(*graph);
		boost::scoped_array<int> vertexDistances(new int[nVertices * nVertices]);
		int* vertexDistancePtr = vertexDistances.get();

		ContextImpl::twoDArray tmp;
		tmp.base = vertexDistances.get();
		tmp.dim = nVertices;

		ContextImpl::constant_property_map_int<Context::inputGraph::edge_descriptor, 1> edgeWeights;
		boost::johnson_all_pairs_shortest_paths(*graph, tmp, boost::weight_map(edgeWeights));
		
		
		edgeDistances = boost::shared_array<int>(new int[nEdges * nEdges]);
		int* edgeDistancePtr = edgeDistances.get();
		memset(edgeDistancePtr, 0, (int)(sizeof(int)*nEdges*nEdges));

		Context::internalGraph::edge_iterator currentFirst, end;
		boost::tie(currentFirst, end) = boost::edges(*graph);
		while(currentFirst != end)
		{
			int firstIndex = boost::get(boost::edge_index, *graph, *currentFirst);

			Context::internalGraph::edge_iterator currentSecond = currentFirst;
			currentSecond++;
			while(currentSecond != end)
			{
				int secondIndex = boost::get(boost::edge_index, *graph, *currentSecond);
				int possibilities[4] = 
				{
					vertexDistancePtr[currentFirst->m_source + nVertices*currentSecond->m_source], 
					vertexDistancePtr[currentFirst->m_source + nVertices*currentSecond->m_target], 
					vertexDistancePtr[currentFirst->m_target + nVertices*currentSecond->m_target], 
					vertexDistancePtr[currentFirst->m_target + nVertices*currentSecond->m_source]
				};
				
				edgeDistancePtr[secondIndex + nEdges * firstIndex] = edgeDistancePtr[firstIndex + nEdges * secondIndex] = *std::min_element(possibilities + 0, possibilities + 4)+1;
				currentSecond++;
			}
			currentFirst++;
		}
	}
	Context::Context()
		:graph(NULL), directedGraph(NULL), vertexPositions(NULL), nEdges(0)
	{}
	Context Context::emptyContext()
	{
		Context result;

		boost::shared_ptr<Context::internalGraph> graph(new Context::internalGraph(2));
		boost::add_edge(0, 1, 0, *graph);
		boost::shared_ptr<const Context::internalGraph> constGraph = boost::static_pointer_cast<const Context::internalGraph>(graph);
		result.graph.swap(constGraph);

		boost::shared_ptr<Context::internalDirectedGraph> directedGraph(new Context::internalDirectedGraph(2));
		boost::add_edge(0, 1, 0, *directedGraph);
		boost::add_edge(1, 0, 1, *directedGraph);
		boost::shared_ptr<const Context::internalDirectedGraph> constDirectedGraph = boost::static_pointer_cast<const Context::internalDirectedGraph>(directedGraph);
		result.directedGraph.swap(constDirectedGraph);

		boost::shared_ptr<std::vector<int> > interestVertices(new std::vector<int>(2));
		(*interestVertices)[0] = 0; (*interestVertices)[1] = 1;
		boost::shared_ptr<const std::vector<int> > constInterestVertices = boost::static_pointer_cast<const std::vector<int> >(interestVertices);
		result.interestVertices.swap(constInterestVertices);

		boost::shared_ptr<std::vector<vertexPosition> > vertexPositions(new std::vector<vertexPosition>(2));
		(*vertexPositions)[0] = vertexPosition(vertexPosition::first_type(0.0), vertexPosition::second_type(0.0)); (*vertexPositions)[1] = vertexPosition(vertexPosition::first_type(10.0), vertexPosition::second_type(0.0));
		boost::shared_ptr<const std::vector<vertexPosition> > constVertexPositions = boost::static_pointer_cast<const std::vector<vertexPosition> >(vertexPositions);
		result.vertexPositions.swap(constVertexPositions);
		
		boost::shared_array<int> edgeDistances(new int[4]);
		edgeDistances[0] = edgeDistances[3] = 0;
		edgeDistances[1] = edgeDistances[2] = 1;
		result.edgeDistances.swap(edgeDistances);

		result.edgeResidualCapacityVector.resize(2);
		result.capacityVector.resize(2);
		result.vertexPredecessorVector.resize(2);
		result.colorVector.resize(2);
		result.distanceVector.resize(2);

		result.nEdges = 0;

		return result;
	}
	Context Context::gridContext(int gridDimension, boost::shared_ptr<const std::vector<int> > interestVertices, capacityDistribution&& distribution, const mpfr_class& threshold)
	{
		if((*interestVertices).size() != 2)
		{
			throw std::runtime_error("Two interest vertices must be entered");
		}

		boost::shared_ptr<std::vector<vertexPosition> > vertexPositions(new std::vector<vertexPosition>(gridDimension * gridDimension));
		boost::shared_ptr<Context::inputGraph> graph(new Context::inputGraph(gridDimension * gridDimension));
		for(int i = 0; i < gridDimension; i++)
		{
			for(int j = 0; j < gridDimension; j++)
			{
				(*vertexPositions)[i + j * gridDimension] = vertexPosition((float)i*100, (float)j*100);
				if(i != gridDimension - 1) boost::add_edge(i + j*gridDimension, i + 1 +j*gridDimension, *graph);
				if(j != gridDimension - 1) boost::add_edge(i + j*gridDimension, i + (j+1)*gridDimension, *graph);
			}
		}
		return Context(graph, interestVertices, vertexPositions, std::move(distribution), threshold);
	}
	Context Context::completeContext(int nVertices, int nInterestVertices, capacityDistribution&& distribution, const mpfr_class& threshold)
	{
		if(nInterestVertices != 2)
		{
			throw std::runtime_error("Input interestVertices must be 2");
		}
		boost::shared_ptr<std::vector<vertexPosition> > vertexPositions(new std::vector<vertexPosition>(nVertices));
		boost::shared_ptr<Context::inputGraph> graph(new Context::inputGraph(nVertices));

		boost::shared_ptr<std::vector<int> > interestVertices(new std::vector<int>(nInterestVertices));
		std::copy(boost::counting_iterator<int>(1), boost::counting_iterator<int>(nInterestVertices), (*interestVertices).begin());

		const double pi = 3.14159265359;
		for(int i = 0; i < nVertices; i++)
		{
			(*vertexPositions)[i] = vertexPosition((float)cos(2*pi*i/nVertices), (float)sin(2*pi*i/nVertices));
			for(int j = i+1; j < nVertices; j++)
			{
				boost::add_edge(i, j, *graph);
			}
		}

		return Context(graph, interestVertices, vertexPositions, std::move(distribution), threshold);
	}
	const Context::internalGraph& Context::getGraph() const
	{
		return *graph;
	}
	const std::vector<int>& Context::getInterestVertices() const
	{
		return *interestVertices;
	}
	const std::vector<Context::vertexPosition>& Context::getVertexPositions() const
	{
		return *vertexPositions;
	}
	Context Context::fromFile(std::string path, bool& successful, boost::shared_ptr<const std::vector<int> > interestVertices, std::string& message, capacityDistribution&& distribution, const mpfr_class& threshold)
	{
		std::ifstream input(path);
		if(!input.is_open())
		{
			successful = false;
			return Context();
		}
		boost::dynamic_properties properties;
		boost::shared_ptr<inputGraph> graph(new inputGraph());

		boost::vector_property_map<int> orderingProperty;
		properties.property("order", orderingProperty);

		boost::vector_property_map<float> xProperty, yProperty;
		properties.property("x", xProperty);
		properties.property("y", yProperty);

		boost::read_graphml(input, *graph, properties);

		boost::shared_ptr<std::vector<vertexPosition> > vertexPositions(new std::vector<vertexPosition>());
		for(auto xIterator = xProperty.storage_begin(), yIterator = yProperty.storage_begin(); xIterator != xProperty.storage_end(); xIterator++, yIterator++)
		{
			vertexPositions->push_back(vertexPosition(*xIterator, *yIterator));
		}

		successful = true;

		int maxInterest = *std::max_element(interestVertices->begin(), interestVertices->end());
		int minInterest = *std::min_element(interestVertices->begin(), interestVertices->end());
		if(minInterest < 0 || maxInterest >= (int)boost::num_vertices(*graph))
		{
			successful = false;
			message = "Invalid vertex indices entered for input interestVertices";
			return Context();
		}
		if(interestVertices->size() != 2)
		{
			throw std::runtime_error("Two interest vertices must be entered");
		}
		return Context(graph, interestVertices, vertexPositions, std::move(distribution), threshold);
	}
	const int* Context::getEdgeDistances() const
	{
		return edgeDistances.get();
	}
	const capacityDistribution& Context::getDistribution() const
	{
		return distribution;
	}
	Context::~Context()
	{
	}
	void Context::constructDirectedGraph()
	{
		boost::shared_ptr<internalDirectedGraph> directedGraph(new internalDirectedGraph(boost::num_vertices(*graph)));
		internalDirectedGraph& directedGraphRef = *directedGraph;
	
		boost::property_map<internalDirectedGraph, boost::edge_reverse_t>::type reverseMap = boost::get(boost::edge_reverse, directedGraphRef);

		internalGraph::edge_iterator start, end, current;
		boost::tie(start, end) = boost::edges(*graph);
		//add edges
		int counter = 0;
		current = start;
		while(current != end)
		{
			std::pair<internalDirectedGraph::edge_descriptor, bool> firstEdgePair = boost::add_edge(current->m_source, current->m_target, counter++, directedGraphRef);

			std::pair<internalDirectedGraph::edge_descriptor, bool> secondEdgePair = boost::add_edge(current->m_target, current->m_source, counter++, directedGraphRef);

			boost::put(reverseMap, firstEdgePair.first, secondEdgePair.first);
			boost::put(reverseMap, secondEdgePair.first, firstEdgePair.first);

			current++;
		}

		boost::shared_ptr<const internalDirectedGraph> constDirectedGraph = boost::static_pointer_cast<const internalDirectedGraph>(directedGraph);
		this->directedGraph.swap(constDirectedGraph);
	}
	const Context::internalDirectedGraph& Context::getDirectedGraph() const
	{
		return *directedGraph;
	}
	std::vector<double>& Context::getCapacityVector() const
	{
		return capacityVector;
	}
	std::size_t Context::getNEdges() const
	{
		return nEdges;
	}
	double Context::getMaxFlow(std::vector<double>& capacities, Context::internalGraph::vertex_descriptor source, Context::internalGraph::vertex_descriptor sink) const
	{
		typedef boost::property_map<Context::internalDirectedGraph, boost::edge_index_t>::const_type edgeIndexMapType;
		typedef boost::property_map<Context::internalDirectedGraph, boost::vertex_index_t>::const_type vertexIndexMapType;
		typedef boost::iterator_property_map<typename std::vector<double>::iterator, edgeIndexMapType> edgeCapacityMapType;
		typedef boost::iterator_property_map<typename std::vector<Context::internalDirectedGraph::edge_descriptor>::iterator, vertexIndexMapType> vertexPredecessorMapType;
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
	double Context::getMaxFlow(std::vector<double>& capacities) const
	{
		return getMaxFlow(capacities, (*interestVertices)[0], (*interestVertices)[1]);
	}
}
