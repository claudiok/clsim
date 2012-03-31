//  tree_suite: Python-land indexing and manipulation for I3Tree
//  (C) 2012 Jakob van Santen <vansanten@wisc.edu>
//           and the IceCube Collaboration
//

#ifndef ICETRAY_PYTHON_TREE_INDEXING_SUITE_HPP_INCLUDED
#define ICETRAY_PYTHON_TREE_INDEXING_SUITE_HPP_INCLUDED

#include <boost/python/iterator.hpp>
#include <boost/python/call_method.hpp>
#include <boost/python/tuple.hpp>

#include <boost/python/object/iterator_core.hpp>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/next_prior.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/erase.hpp>
#include <boost/mpl/replace.hpp>

namespace bp = boost::python;

namespace boost { namespace python {

namespace tree_suite_detail {
	template <typename Sig>
	struct LastArg
	{
		typedef typename mpl::deref<typename mpl::prior<
		    typename mpl::end<Sig>::type>::type>::type type;
	};
	template <typename Sig>
	struct ReturnType
	{
		typedef typename mpl::deref<
		    typename mpl::begin<Sig>::type>::type type;
	};
}

template <class Container>
class tree_indexing_suite : public bp::def_visitor<tree_indexing_suite<Container > >{
public:
	
	typedef typename Container::value_type value_type;
	typedef typename Container::pre_order_iterator pre_order_iterator;
	typedef typename Container::post_order_iterator post_order_iterator;
	typedef typename Container::sibling_iterator sibling_iterator;
		
	static void
	index_translation_error(int offset, unsigned depth) throw(python::error_already_set)
	{
		std::ostringstream os;
		os << "Index (" << offset << ") out of range at level " << depth;
		PyErr_SetString(PyExc_IndexError, os.str().c_str());
		throw python::error_already_set();
	}
	
	static typename Container::sibling_iterator
	translate_index(Container &self, bp::tuple &index)
	{
		unsigned maxdepth = len(index)-1;
		sibling_iterator it = self.begin();
		for (unsigned depth=0; depth < maxdepth; depth++) {
			// Advance among the siblings by the given index
			int offset = extract<int>(index[depth]);
			it += offset;
			
			if (!self.is_valid(it))
				index_translation_error(offset, depth);
			// Drop down a level at the indicated position
			it = self.begin(it);
		}
		
		int offset = extract<int>(index[maxdepth]);
		it += offset;
		
		if (!self.is_valid(it))
			index_translation_error(offset, maxdepth);
		return it;
	}
	
	static bp::tuple
	translate_iterator(Container &self, sibling_iterator it)
	{
		std::vector<int> indices(self.depth(it)+1, 0);
		for (unsigned depth = indices.size()-1; depth > 0; depth--) {
			indices[depth] = self.index(it);
			it = self.parent(it);
		}
		
		indices[0] = self.index(it);
			
		return bp::tuple(indices);
	}
	
	// A basic range wrapper much like boost::python::range
	template <typename Iterator>
	class iterator : public Iterator {
	public:	
		iterator(Container &target, Iterator &begin, Iterator &end, bool addr)
		    : self_(target), begin_(begin), end_(end), return_address_(addr) {};
		python::object next()
		{
			if (begin_ == end_)
				objects::stop_iteration_error();
			if (return_address_) {
				python::tuple addr = translate_iterator(self_, begin_);
				return python::make_tuple(addr, *begin_++);
			} else
				return python::object(*begin_++);
		}		
	private:
		Container &self_;
		Iterator begin_, end_;
		bool return_address_;
	};
	
	template <typename Iterator>
	static void
	register_iterator(const char *name)
	{
		typedef iterator<Iterator> IteratorWrapper;
		class_<IteratorWrapper>(name, bp::no_init)
			.def("__iter__", objects::identity_function())
			.def("next", &IteratorWrapper::next)
		;
	}
	
	template <typename Accessor>
	static inline iterator<Accessor>
	make_iterator(Container &x, Accessor begin, Accessor end, bool addr)
	{
		return iterator<Accessor>(x, begin, end, addr);
	}
	
	static iterator<pre_order_iterator>
	get_iterator(Container &x)
	{
		return make_iterator(x, x.begin(), x.end(), false);
	};
	
	static iterator<pre_order_iterator>
	get_keyed_iterator(Container &x)
	{
		return make_iterator(x, x.begin(), x.end(), true);
	};
	
	static iterator<post_order_iterator>
	get_post_order_iterator(Container &x, bool addr)
	{
		return make_iterator(x, x.begin_post(), x.end_post(), addr);
	};
	
	static iterator<sibling_iterator>
	get_sibling_iterator(Container &x, tuple pos, bool addr)
	{
		typename Container::sibling_iterator it = translate_index(x, pos);
		return make_iterator(x, x.begin(it), x.end(it), addr);
	};
	
	static bp::tuple
	add_root(Container &self, const value_type &v)
	{
		if (self.size() == 0)
			return translate_iterator(self, self.set_head(v));
		else
			return translate_iterator(self, self.insert(self.begin(), v));
	}
	
	static value_type
	getitem(Container &self, bp::tuple index)
	{
		return *translate_index(self, index);
	}
	
	static value_type
	getitem_simple(Container &self, int index)
	{
		return getitem(self, make_tuple(index));
	}
	
	static void
	setitem(Container &self, bp::tuple index, const value_type &v)
	{
		*translate_index(self, index) = v;
	}
	
	static void
	setitem_simple(Container &self, int index, const value_type &v)
	{
		setitem(self, make_tuple(index), v);
	}
	
	static bp::list
	argfind(Container &self, bp::object comparator)
	{
		bp::list results;
		typename Container::iterator it;
		for (it = self.begin(); it != self.end(); it++)
			if (bp::extract<bool>(comparator(*it)))
				results.append(translate_iterator(self, it));
		return results;
	}
	
	// iterator_adapter<F,S>::make_function(F f) takes a member-function pointer
	// that takes and returns tree iterators, and automatically generates a wrapper
	// function that converts an address tuple to the equivalent iterator (e.g. (0,1,0)
	// to first_root->second_child->first_child), calls the function, and converts
	// the result back to an address tuple.  
	template <typename F, typename S>
	struct iterator_adapter : public bp::object {
		
		typedef F FunctionPtr;
		typedef S HeldSignature;
		typedef typename tree_suite_detail::ReturnType<S>::type ReturnType;
		typedef typename tree_suite_detail::LastArg<S>::type Arg1;
		typedef typename mpl::replace<S, ReturnType, bp::tuple>::type WrapperSignature;		
		
		typedef iterator_adapter<F,S> Adapter;
		
		FunctionPtr f_;
		
		iterator_adapter<F,S>(F f) : f_(f) {};
		
		// SFINAE dispatch to specializations for different signatures
		template <typename S1>
		static typename bp::object
		make_function(F f, S1 s)
		{
			return Adapter::make_function_impl(f, s);
		}
		
	private:
		
		// Methods of the form: iterator (Container::*)(iterator, value_type)
		template <typename S1>
		static typename enable_if<
		    mpl::equal_to<typename mpl::size<S1>::type,
			          typename mpl::long_<4l>::type>,
		    bp::object>::type
		make_function_impl(F f, S1 s)
		{
			return bp::make_function(boost::bind(&Adapter::call_inserter, Adapter(f), _1, _2, _3),
			    default_call_policies(), args("self", "index", "value"), WrapperSignature());
		}
		
		
		// Methods of the form: iterator (Container::*)(iterator)
		template <typename S1>
		static typename enable_if<typename mpl::and_<
			                  typename mpl::equal_to<typename mpl::size<S1>::type,
					                         typename mpl::long_<3l>::type>::type,
				          typename is_same<typename mpl::at<S1, typename mpl::long_<2> >::type,
					                   ReturnType>::type>,
		    bp::object>::type
		make_function_impl(F f, S1 s)
		{
			return bp::make_function(boost::bind(&Adapter::call_single_iter, Adapter(f), _1, _2),
			    default_call_policies(), args("self", "index"), WrapperSignature());
			
		}
		
		bp::tuple call_inserter(const Container &self, bp::tuple index, Arg1 a1)
		{
			return translate_call(self, boost::bind(f_, const_cast<Container*>(&self), _1, a1), index);
		}
		
		bp::tuple call_single_iter(const Container &self, bp::tuple index)
		{
			return translate_call(self, boost::bind(f_, const_cast<Container*>(&self), _1), index);
		}
		
		// Translate tuple to iterator, call functor, translate result back
		template <typename Caller>
		static inline bp::tuple
		translate_call(const Container &selfc, Caller caller, bp::tuple index)
		{
			Container &self = const_cast<Container&>(selfc);
			ReturnType it = translate_index(self, index);
			it = caller(it);
			if (self.is_valid(it))
				return translate_iterator(self, it);
			else
				return bp::make_tuple(0);
		}
		
	};
	
	template <typename F>
	static bp::object
	make_adapter(F f)
	{
		return make_adapter_impl(f, detail::get_signature(f));
	}
	
	template <typename F, typename Sig>
	static bp::object
	make_adapter_impl(F f, Sig sig)
	{
		return iterator_adapter<F,Sig>::make_function(f, sig);
	}
	
	template <class Class>
	static void
	visit(Class& cl)
	{
		typedef typename mpl::if_<
		    is_class<value_type>
		    , return_internal_reference<>
		    , default_call_policies
		    >::type get_data_return_policy;
		
		bp::scope outer = cl;
		
		register_iterator<pre_order_iterator>("pre_order_iterator");
		register_iterator<post_order_iterator>("post_order_iterator");
		register_iterator<sibling_iterator>("sibling_iterator");
		
		typedef sibling_iterator (Container::*inserter_func)(sibling_iterator, const value_type&);
		typedef sibling_iterator (Container::*itermanip_func)(sibling_iterator);		
				
		cl
			.def("__len__", &Container::size)
			.def("__iter__", &get_iterator)
			.def("__getitem__", &getitem)
			.def("__getitem__", &getitem_simple)
			.def("__setitem__", &setitem)
			.def("__setitem__", &setitem_simple)
			.def("iteritems", &get_keyed_iterator)
			.def("iter_children", &get_sibling_iterator,     arg("self"), arg("index"), arg("return_indices")=true)
			.def("iter_postorder", &get_post_order_iterator, arg("self"), arg("return_indices")=true)
			.def("argfind", &argfind, arg("self"), arg("keyfunc"))
			.def("add_root", &add_root, arg("self"), arg("value"))
			.def("insert",       make_adapter((inserter_func)&Container::insert))
			.def("insert_after", make_adapter((inserter_func)&Container::insert_after))
			.def("replace",      make_adapter((inserter_func)&Container::replace))
			.def("append_child", make_adapter((inserter_func)&Container::append_child))
			.def("erase",        make_adapter((itermanip_func)&Container::erase))
		;
		
	}
	
};


} } // namespace boost::python

#endif
