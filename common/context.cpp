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
		source = other.source;
		sink = other.sink;
		distributions = std::move(other.distributions);

		edgeResidualCapacityVector.swap(other.edgeResidualCapacityVector);
		capacityVector.swap(other.capacityVector);
		directedGraph.swap(other.directedGraph);
		vertexPredecessorVector.swap(other.vertexPredecessorVector);
		colorVector.swap(other.colorVector);
		distanceVector.swap(other.distanceVector);

		std::swap(threshold, other.threshold);
		return *this;
	}
	Context::Context(Context&& other)
	{
		graph.swap(other.graph);
		source = other.source;
		sink = other.sink;
		distributions = std::move(other.distributions);

		edgeResidualCapacityVector.swap(other.edgeResidualCapacityVector);
		capacityVector.swap(other.capacityVector);
		directedGraph.swap(other.directedGraph);
		vertexPredecessorVector.swap(other.vertexPredecessorVector);
		colorVector.swap(other.colorVector);
		distanceVector.swap(other.distanceVector);

		std::swap(threshold, other.threshold);
	}
	Context::Context(boost::shared_ptr<const inputGraph> unorderedGraph, int source, int sink, std::vector<capacityDistribution>&& distributions, const mpfr_class& threshold)
		: source(source), sink(sink), distributions(std::move(distributions))
	{
		std::size_t nVertices = boost::num_vertices(*unorderedGraph);
		int minVertexIndex = std::min(source, sink);
		int maxVertexIndex = std::max(source, sink);
		if(minVertexIndex < 0 || maxVertexIndex >= (int)nVertices)
		{
			throw std::runtime_error("Input interestVertices was out of range");
		}
		std::size_t nEdges = boost::num_edges(*unorderedGraph);
		if(this->distributions.size() != nEdges)
		{
			throw std::runtime_error("Input distributions must have an entry for every edge");
		}
		
		//Double the number of edges because these are going to be used on the directed version
		edgeResidualCapacityVector.resize(2*nEdges);
		vertexPredecessorVector.resize(2*nEdges);
		capacityVector.resize(2*nEdges);
		vertexPredecessorVector.resize(2*nEdges);
		colorVector.resize(2*nEdges);
		distanceVector.resize(2*nEdges);

		if(boost::connected_components(*unorderedGraph, &(capacityVector[0])) > 1)
		{
			throw std::runtime_error("Input graph must contain a single connected component");
		}

		//Construct a new graph, where the edge indices are assigned consecutively
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
	}
	Context::Context()
		:graph(NULL), directedGraph(NULL)
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

		result.source = 0;
		result.sink = 1;
		
		result.edgeResidualCapacityVector.resize(2);
		result.capacityVector.resize(2);
		result.vertexPredecessorVector.resize(2);
		result.colorVector.resize(2);
		result.distanceVector.resize(2);

		return result;
	}
	Context Context::gridContext(int gridDimension, int source, int sink, std::vector<capacityDistribution>&& distributions, const mpfr_class& threshold)
	{
		boost::shared_ptr<Context::inputGraph> graph(new Context::inputGraph(gridDimension * gridDimension));
		for(int i = 0; i < gridDimension; i++)
		{
			for(int j = 0; j < gridDimension; j++)
			{
				if(i != gridDimension - 1) boost::add_edge(i + j*gridDimension, i + 1 +j*gridDimension, *graph);
				if(j != gridDimension - 1) boost::add_edge(i + j*gridDimension, i + (j+1)*gridDimension, *graph);
			}
		}
		return Context(graph, source, sink, std::move(distributions), threshold);
	}
	Context Context::completeContext(int nVertices, std::vector<capacityDistribution>&& distributions, const mpfr_class& threshold)
	{
		boost::shared_ptr<Context::inputGraph> graph(new Context::inputGraph(nVertices));

		for(int i = 0; i < nVertices; i++)
		{
			for(int j = i+1; j < nVertices; j++)
			{
				boost::add_edge(i, j, *graph);
			}
		}

		return Context(graph, 0, 1, std::move(distributions), threshold);
	}
	const Context::internalGraph& Context::getGraph() const
	{
		return *graph;
	}
	int Context::getSource() const
	{
		return source;
	}
	int Context::getSink() const
	{
		return sink;
	}
	Context Context::fromFile(std::string path, bool& successful, int source, int sink, std::string& message, std::vector<capacityDistribution>&& distributions, const mpfr_class& threshold)
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

		successful = true;

		if(std::min(source, sink) < 0 || std::max(source, sink) >= (int)boost::num_vertices(*graph))
		{
			successful = false;
			message = "Invalid vertex indices entered for input interestVertices";
			return Context();
		}
		return Context(graph, source, sink, std::move(distributions), threshold);
	}
	const capacityDistribution& Context::getDistribution(int edgeIndex) const
	{
		return distributions[edgeIndex];
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
/*	double Context::getMaxFlow(std::vector<double>& capacities, Context::internalGraph::vertex_descriptor source, Context::internalGraph::vertex_descriptor sink) const
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
	}*/
}
