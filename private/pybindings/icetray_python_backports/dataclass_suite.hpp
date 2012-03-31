//  dataclass_suite: provide all the goodies dataclasses should have in
//                   a single visitor, including:
//  - indexing (for maps, vectors, and trees)
//  - copy constructors
//  - __str__ representations
//  - boost::serialization-backed pickling
//  - disabling __setattr__ to avoid API confusion
//  
//  (C) 2012 Jakob van Santen <vansanten@wisc.edu>
//           and the IceCube Collaboration



#include <boost/python/def_visitor.hpp>

#include "copy_suite.hpp"
#include "stream_to_string.hpp"
#include "boost_serializable_pickle_suite.hpp"
#include "std_map_indexing_suite.hpp"
#include "std_vector_indexing_suite.hpp"
#include "tree_indexing_suite.hpp"

#include <boost/utility/enable_if.hpp>

#include <boost/mpl/if.hpp>
#include <boost/mpl/or.hpp>

#include <boost/static_assert.hpp>
#include <boost/typeof/typeof.hpp>

namespace boost { namespace python {

namespace detail {

#define HAS_TYPEDEF(def, name)                                           \
	template<typename T>                                             \
	struct name                                                      \
	{                                                                \
	private:                                                         \
	    typedef char                      yes;                       \
	    typedef struct { char array[2]; } no;                        \
                                                                         \
	    template<typename C> static yes test(typename C::def*);      \
	    template<typename C> static no  test(...);                   \
	public:                                                          \
	    static const bool value = sizeof(test<T>(0)) == sizeof(yes); \
	}                                                                \

HAS_TYPEDEF(iterator, has_iterator);
HAS_TYPEDEF(fixed_depth_iterator, is_tree);
HAS_TYPEDEF(key_type, is_map);

#undef HAS_TYPEDEF

template <bool B, typename T>
struct has_random_access_iterator_impl;

template<typename T>
struct has_random_access_iterator_impl<false, T>
{
	static bool const value = false;
};

template<typename T>
struct has_random_access_iterator_impl<true, T>
{
	static const bool value = boost::is_same<typename std::iterator_traits<typename T::iterator>::iterator_category,
		    std::random_access_iterator_tag>::value;

};

template <typename T>
struct is_vector {
	static const bool value = has_random_access_iterator_impl<has_iterator<T>::value, T>::value;
};

namespace has_operator {
    typedef char yes;
    typedef struct { char array[2]; } no;
    
    struct anyx { template <class T> anyx(const T &); };
    no operator << (const anyx &, const anyx &);
    
    template <class T> yes check(T const&);
    no check(no);
    
    template <typename StreamType, typename T>
    struct insertion {
        static StreamType & stream;
        static T & x;
        static const bool value = sizeof(check(stream << x)) == sizeof(yes);
    };
}

}
	

	
template <class T>
class dataclass_suite : public bp::def_visitor<dataclass_suite<T > > {
private:
	
	template <class Class, typename U>
	static
	typename boost::enable_if_c<detail::is_map<U>::value>::type
	add_indexing(Class &cl)
	{
		cl.def(std_map_indexing_suite<U>());
	}

	template <class Class, typename U>
	static
	typename boost::enable_if_c<detail::is_vector<U>::value>::type
	add_indexing(Class &cl)
	{
		cl.def(std_vector_indexing_suite<U>());
	}
	
	template <class Class, typename U>
	static
	typename boost::enable_if_c<detail::is_tree<U>::value>::type
	add_indexing(Class &cl)
	{
		cl.def(tree_indexing_suite<U>());
	}

	template <class Class, typename U>
	static
	typename boost::disable_if<boost::mpl::or_<
		boost::mpl::bool_<detail::is_map<U>::value>,
		boost::mpl::bool_<detail::is_vector<U>::value>,
		boost::mpl::bool_<detail::is_tree<U>::value> > >::type
	add_indexing(Class &cl) {}
	
	template <class Class, typename U>
	static
	typename boost::enable_if_c<detail::has_operator::insertion<std::ostream, U>::value>::type
	add_string_to_stream(Class &cl)
	{
		cl.def("__str__", &stream_to_string<U>);
	}

	template <class Class, typename U>
	static
	typename boost::disable_if_c<detail::has_operator::insertion<std::ostream, U>::value>::type
	add_string_to_stream(Class &cl) {}
	
public:
	template <class Class>
	static void
	visit(Class& cl)
	{
		// Add any must-have visitors here. Any that rely on non-universal
		// class attributes (like the indexing suites) should be
		// conditionalized to fail gracefully, following the SFINAE
		// pattern above.
		cl.def(copy_suite<T>());
		cl.def_pickle(boost_serializable_pickle_suite<T>());
		add_string_to_stream<Class, T>(cl);
		add_indexing<Class, T>(cl);
		cl.def(freeze());
	}
	
};

}}