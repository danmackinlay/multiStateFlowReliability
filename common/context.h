#ifndef CONTEXT_HEADER_GUARD
#define CONTEXT_HEADER_GUARD
#include <boost/graph/adjacency_list.hpp>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <vector>
#include "includeMPFR.h"
#include "capacityDistribution.h"
namespace multistateTurnip
{
	class Context : public boost::noncopyable
	{
	public:
		friend class boost::serialization::access;
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> inputGraph;
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, boost::property<boost::edge_index_t, int> > internalGraph;
		typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::bidirectionalS> directedGraphTraits;
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property, 
			boost::property<boost::edge_index_t, int, 
				boost::property<boost::edge_reverse_t, directedGraphTraits::edge_descriptor>
			> 
		> internalDirectedGraph;

		Context(boost::shared_ptr<const inputGraph> graph, int source, int sink, std::vector<capacityDistribution>&& distributions, const mpfr_class& threshold);

		Context& operator=(Context&& other);
		Context(Context&& other);
		~Context();
		
		static Context gridContext(int gridDimension, int source, int sink, std::vector<capacityDistribution>&& distributions, const mpfr_class& threshold);
		static Context fromFile(std::string path, bool& successful, int source, int sink, std::string& message, std::vector<capacityDistribution>&& distributions, const mpfr_class& threshold);
		static Context emptyContext();
		const internalGraph& getGraph() const;
		const internalDirectedGraph& getDirectedGraph() const;
		int getSource() const;
		int getSink() const;
		const capacityDistribution& getDistribution(int edgeIndex) const;
		std::vector<double>& getCapacityVector() const;
		static Context completeContext(int nVertices, std::vector<capacityDistribution>&& distributions, const mpfr_class& threshold);
//		double getMaxFlow(std::vector<double>& capacities) const;
//		double getMaxFlow(std::vector<double>& capacities, Context::internalGraph::vertex_descriptor source, Context::internalGraph::vertex_descriptor sink) const;
	private:
		Context& operator=(Context const& other);
		Context();
		void constructDirectedGraph();
		boost::shared_ptr<const internalGraph> graph;
		boost::shared_ptr<const internalDirectedGraph> directedGraph;
		int source, sink;

		std::vector<capacityDistribution> distributions;
		mpfr_class threshold;
		
		//These are used in the min paths call. Stored here so they can be reused. 
		mutable std::vector<double> edgeResidualCapacityVector;
		mutable std::vector<double> capacityVector;
		mutable std::vector<internalDirectedGraph::edge_descriptor> vertexPredecessorVector;
		mutable std::vector<boost::default_color_type> colorVector;
		mutable std::vector<int> distanceVector;
	};
}
namespace boost
{
	namespace networkReliabilityImpl
	{
		struct customStruct
		{};
	}
	//Fix screw ups in boykov_kolmogorov_max_flow (Note that it doesn't take the edge_capacity_t from the named parameters, only from the graph).
	//May be fixed by boost at some point, in which case delete this. 
	template<> struct property_map<const multistateTurnip::Context::internalDirectedGraph, edge_capacity_t>
	{
		typedef networkReliabilityImpl::customStruct const_type;
	};
	template<> struct property_traits<networkReliabilityImpl::customStruct>
	{
		typedef double value_type;
	};
}
#endif
