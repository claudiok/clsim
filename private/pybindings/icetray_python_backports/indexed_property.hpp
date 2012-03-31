//  indexed_property: expose getter/setter functions with multiple 
//                    arguments as proxy objects with __getitem__/__setitem__.
//  
//  (C) 2012 Jakob van Santen <vansanten@wisc.edu>
//           and the IceCube Collaboration

#ifndef ICETRAY_PYTHON_INDEXED_PROPERTY_HPP_INCLUDED
#define ICETRAY_PYTHON_INDEXED_PROPERTY_HPP_INCLUDED

#include <boost/python/def_visitor.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/signature.hpp>

#include <boost/mpl/at.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/pop_front.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/greater.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/advance.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>

namespace boost { namespace python {

template <typename G, typename GS, typename S, typename SS>
class indexed_property_impl {

private:
	typedef G Getter;
	typedef S Setter;
	
	// Getter must have at least one argument, else we wouldn't need to jump through this
	// many hoops to get a nice property interface. The current implementation supports up
	// to 3 key parameters.
	BOOST_MPL_ASSERT(( mpl::greater<typename mpl::size<GS>::type, typename mpl::long_<2l>::type > ));
	BOOST_MPL_ASSERT(( mpl::less<typename mpl::size<GS>::type, typename mpl::long_<6l>::type > ));
		
	typedef typename mpl::at<GS, typename mpl::long_<0>::type>::type ValueType;
	typedef typename mpl::at<GS, typename mpl::long_<1>::type>::type Target;
	typedef typename mpl::pop_front<typename mpl::pop_front<GS>::type>::type GetArgs;
		
	// Test if we actually have a setter
	typedef typename mpl::not_<typename mpl::bool_<is_same<void (*)(), S>::value > > have_setter;	
	typedef typename mpl::if_<have_setter, 
		typename mpl::push_back<GetArgs, ValueType>::type, void*>::type DeducedSetArgs;
	typedef typename mpl::if_<have_setter, 
		typename mpl::pop_front<typename mpl::pop_front<SS>::type>::type, void*>::type SetArgs;
	typedef typename mpl::if_<have_setter, 
		typename mpl::at<SS, typename mpl::long_<0>::type>::type, void>::type SetReturnType;
	
	// Ensure that the Setter is of the form void (Target::*)(..., ValueType)
	// FIXME: how do I compare synthesized MPL sequence types?
	// BOOST_MPL_ASSERT(( is_same<DeducedSetArgs, SetArgs> ));
	BOOST_MPL_ASSERT(( is_same<SetReturnType, void> ));
	
	Getter getter_;
	Setter setter_;
	
	// Convert reference types to value types
	template <typename T> struct value_type { typedef T type; };
	template <typename T> struct value_type<T&> { typedef T type; };
	typedef shared_ptr<typename value_type<Target>::type> TargetPtr;
	
	// Make getters to return references if needed
	typedef typename mpl::and_<
		is_reference<ValueType>,
		is_class<typename value_type<ValueType>::type>
		>::type valuetype_is_reference;
	
	typedef typename mpl::if_<valuetype_is_reference,
	    return_internal_reference<1>,
		default_call_policies
	    >::type getter_policy;	
	
	struct name_collector {
		name_collector(std::string &n) : names(n) {};
		std::string &names;
		template <typename T>
		void operator()(T t)
		{
			if (names.size() > 0)
				names += std::string(", ") + type_id<T>().name();
			else
				names += type_id<T>().name();
		}
	};
	
	template <typename Args>
	static std::string
	stringify_signature()
	{
		std::string names;
		name_collector coll(names);
		mpl::for_each<Args>(coll);
		return names;
	}
	
	template <typename Args, unsigned N>
	struct invoke_getter;
	
	template <typename Args>
	struct invoke_getter<Args,2> {
		typedef typename mpl::at<Args, mpl::long_<0l> >::type A0;
		typedef typename mpl::at<Args, mpl::long_<1l> >::type A1;
		
		static ValueType apply(TargetPtr target, G g, python::tuple args)
		{
			return (*target.*g)(
			    python::extract<A0>(args[0]),
			    python::extract<A1>(args[1])
				);
		}
	};
	
	template <typename Args>
	struct invoke_getter<Args,3> {
		typedef typename mpl::at<Args, mpl::long_<0l> >::type A0;
		typedef typename mpl::at<Args, mpl::long_<1l> >::type A1;
		typedef typename mpl::at<Args, mpl::long_<2l> >::type A2;
		
		static ValueType apply(TargetPtr target, G g, python::tuple args)
		{
			
			A0 a0 = python::extract<A0>(args[0]);
			A1 a1 = python::extract<A1>(args[1]);
			A2 a2 = python::extract<A1>(args[2]);
		
			return (*target.*g)(a0, a1, a2);
		}
	};
	
	template <typename Args, unsigned N>
	struct invoke_setter;
	
	template <typename Args>
	struct invoke_setter<Args,3> {
		typedef typename mpl::at<Args, mpl::long_<0l> >::type A0;
		typedef typename mpl::at<Args, mpl::long_<1l> >::type A1;
		
		static void apply(TargetPtr target, S s, python::tuple args, ValueType v)
		{
			return (*target.*s)(
			    python::extract<A0>(args[0]),
			    python::extract<A1>(args[1]),
				v
				);
		}
	};
	
	template <typename Args>
	struct invoke_setter<Args,4> {
		typedef typename mpl::at<Args, mpl::long_<0l> >::type A0;
		typedef typename mpl::at<Args, mpl::long_<1l> >::type A1;
		typedef typename mpl::at<Args, mpl::long_<2l> >::type A2;
		
		static void apply(TargetPtr target, S s, python::tuple args, ValueType v)
		{
			return (*target.*s)(
			    python::extract<A0>(args[0]),
			    python::extract<A1>(args[1]),
				python::extract<A2>(args[2]),	
				v
				);
		}
	};
	
	struct Proxy {
		TargetPtr held_;
		Getter getter_;
		Setter setter_;
		
		Proxy(TargetPtr target, Getter g) : held_(target), getter_(g) {}
		Proxy(TargetPtr target, Getter g, Setter s) : held_(target), getter_(g), setter_(s) {}
		
		static Proxy
		create(TargetPtr target, Getter g)
		{
			return Proxy(target, g);
		}
		
		static Proxy
		create(TargetPtr target, Getter g, Setter s)
		{
			return Proxy(target, g, s);
		}
			
		static object
		make_function(Getter g)
		{
			return python::make_function(bind(&Proxy::create, _1, g), default_call_policies(),
			    mpl::vector<Proxy, TargetPtr>());
		}
		
		static object
		make_function(Getter g, Setter s)
		{
			return python::make_function(bind(&Proxy::create, _1, g, s), default_call_policies(),
			    mpl::vector<Proxy, TargetPtr>());
		}

		// Specializations for arity
		template <typename Args>
		typename enable_if<
		    mpl::equal_to<typename mpl::size<Args>::type,
			          typename mpl::long_<1l>::type>,
		    ValueType>::type
		getitem(typename mpl::at<Args, mpl::long_<0l> >::type a1)
		{
			return (*held_.*getter_)(a1);
		}
		
		template <typename Args>
		typename enable_if<
		    mpl::greater<typename mpl::size<Args>::type,
			          typename mpl::long_<1l>::type>,
		    ValueType>::type
		getitem(python::tuple args)
		{
			return invoke_getter<Args, mpl::size<Args>::type::value>::apply(held_, getter_, args);
		}
		
		template <typename Args>
		typename enable_if<
		    mpl::equal_to<typename mpl::size<Args>::type,
			          typename mpl::long_<2l>::type>,
		    void>::type
  		setitem(typename mpl::at<Args, mpl::long_<0l> >::type a1,
  		        ValueType v)
		{
			return (*held_.*setter_)(a1, v);
		}
		
		template <typename Args>
		typename enable_if<
		    mpl::greater<typename mpl::size<Args>::type,
			          typename mpl::long_<2l>::type>,
		    void>::type
  		setitem(python::tuple args, ValueType v)
		{
			return invoke_setter<Args, mpl::size<Args>::type::value>::apply(held_, setter_, args, v);			
		}
		
	};
	
	struct ReadOnly {};
	struct ReadWrite {};

	template <class Kind, class Class>
	void
	bind_property(Class &cl, class_<Proxy> &proxy_class, const char* propname, Kind k) const;
	
	template <class Class>
	void
	bind_property(Class &cl, class_<Proxy> &proxy_class, const char* propname, ReadOnly k) const
	{
		cl.add_property(propname, Proxy::make_function(getter_));
	}
	
	template <class Class>
	void
	bind_property(Class &cl, class_<Proxy> &proxy_class, const char* propname, ReadWrite k) const
	{
		cl.add_property(propname, Proxy::make_function(getter_, setter_));
		std::string usage = std::string("Usage: ") + std::string(propname) + "["
		    + stringify_signature<SetArgs>() + "] = " + type_id<ValueType>().name();
		proxy_class.def("__setitem__", &Proxy::template setitem<SetArgs>, usage.c_str());
	}
	
public:
	indexed_property_impl<G,GS,S,SS>(G g, S s) : getter_(g), setter_(s) {};
	
	template <class Class, class Options>
	void visit(Class& cl, const char *name, Options &options) const
	{
		
		scope inner = cl;
		
		std::string proxy_name = "_" + std::string(name) + "_proxy";
				
		class_<Proxy> proxy_class(proxy_name.c_str(), no_init);
		proxy_class.def_readonly("parent", &Proxy::held_);
		
		std::string usage = std::string("Usage: ") + name + "[" + stringify_signature<GetArgs>() + "]";
		
		proxy_class.def("__getitem__", &Proxy::template getitem<GetArgs>, usage.c_str(), getter_policy());

		typedef is_same<void (*)(), Setter> no_setter;
		typedef typename mpl::if_c<no_setter::value, ReadOnly, ReadWrite>::type property_type;
		bind_property(cl, proxy_class, name, property_type());
		
	}

};

template <class Class, class Options, typename G, typename GS, typename S, typename SS>
void dispatch_visit(Class& cl, const char *name, Options &options, G g, GS gs, S s, SS ss)
{
	indexed_property_impl<G,GS,S,SS>(g,s).visit(cl, name, options);
}

template <typename G, typename S>
class indexed_property : public def_visitor<indexed_property<G,S> > {
 
public:
	
	indexed_property(G g, S s) : g_(g), s_(s) {};
	
	template <class Class, class Options>
	void visit(Class& cl, const char *name, Options &options) const
	{
		// Kick down to an implementation that knows the function signature
		dispatch_visit(cl, name, options, g_, detail::get_signature(g_), s_, detail::get_signature(s_));
	}
	
private:
	G g_;
	S s_;	
	
};

// Dispatcher for read-only properties
template <typename G>
indexed_property<G, void (*)()>
make_indexed_property(G g)
{
	typedef void (*nullary)();
	nullary nil;
	return indexed_property<G, nullary>(g, nil);
}

// Dispatcher for read-write properties
template <typename G, typename S>
indexed_property<G,S>
make_indexed_property(G g, S s)
{
	return indexed_property<G,S>(g, s);
}
	
}} // namespace boost::python

#endif

