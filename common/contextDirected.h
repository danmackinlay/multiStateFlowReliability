#ifndef CONTEXT_DIRECTED_HEADER_GUARD
#define CONTEXT_DIRECTED_HEADER_GUARD
#include <boost/graph/adjacency_list.hpp>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <vector>
#include "includeMPFR.h"
#include "capacityDistribution.h"
#include "context.h"
namespace multistateTurnip
{
	class ContextDirected : public boost::noncopyable
	{
	public:
		friend class boost::serialization::access;
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, boost::no_property, boost::property<boost::edge_index_t, int> > inputGraph;
		typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS> directedGraphTraits;
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, boost::no_property, boost::property<boost::edge_index_t, int, boost::property<boost::edge_reverse_t, directedGraphTraits::edge_descriptor> > > internalGraph;

		ContextDirected(boost::shared_ptr<const inputGraph> graph, int source, int sink, const std::vector<capacityDistribution>& distributions, const mpfr_class& threshold);

		ContextDirected& operator=(ContextDirected&& other);
		ContextDirected(ContextDirected&& other);
		~ContextDirected();
		
		static ContextDirected gridContextDirected(int gridDimension, int source, int sink, const std::vector<capacityDistribution>& distributions, const mpfr_class& threshold);
		static ContextDirected fromFile(std::string path, bool& successful, int source, int sink, std::string& message, const std::vector<capacityDistribution>& distributions, const mpfr_class& threshold);
		static ContextDirected emptyContextDirected();
		const internalGraph& getGraph() const;
		int getSource() const;
		int getSink() const;
		const capacityDistribution& getDistribution(int edgeIndex) const;
		std::vector<double>& getCapacityVector() const;
		static ContextDirected completeContextDirected(int nVertices, const std::vector<capacityDistribution>& distributions, const mpfr_class& threshold);
//		double getMaxFlow(std::vector<double>& capacities) const;
//		double getMaxFlow(std::vector<double>& capacities, ContextDirected::internalGraph::vertex_descriptor source, ContextDirected::internalGraph::vertex_descriptor sink) const;
	private:
		ContextDirected& operator=(ContextDirected const& other);
		ContextDirected();
		boost::shared_ptr<const internalGraph> graph;
		int source, sink;

		std::vector<capacityDistribution> distributions;
		mpfr_class threshold;
		
		//These are used in the min paths call. Stored here so they can be reused. 
		mutable std::vector<double> edgeResidualCapacityVector;
		mutable std::vector<double> capacityVector;
		mutable std::vector<internalGraph::edge_descriptor> vertexPredecessorVector;
		mutable std::vector<boost::default_color_type> colorVector;
		mutable std::vector<int> distanceVector;
	};
}
namespace boost
{
	//Fix screw ups in boykov_kolmogorov_max_flow (Note that it doesn't take the edge_capacity_t from the named parameters, only from the graph).
	//May be fixed by boost at some point, in which case delete this. 
	template<> struct property_map<const multistateTurnip::ContextDirected::internalGraph, edge_capacity_t>
	{
		typedef networkReliabilityImpl::customStruct const_type;
	};
}
#endif
