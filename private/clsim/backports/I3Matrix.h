#ifndef I3CLSIM_I3MATRIX_H_INCLUDED
#define I3CLSIM_I3MATRIX_H_INCLUDED

#include <icetray/I3FrameObject.h>
#include <boost/numeric/ublas/matrix.hpp>

class I3Matrix : public I3FrameObject,
    public boost::numeric::ublas::matrix<double> {
public:
	typedef boost::numeric::ublas::matrix<double> base;
	typedef base::size_type size_type;
	typedef base::difference_type difference_type;
	typedef base::value_type value_type;
	typedef base::const_reference const_referece;
	typedef base::reference reference;
	typedef base::array_type array_type;
	typedef base::const_closure_type const_closure_type;
	typedef base::closure_type closure_type;
	typedef base::vector_temporary_type vector_temporary_type;
	typedef base::matrix_temporary_type matrix_temporary_type;
	typedef base::storage_category storage_category;
	typedef base::orientation_category orientation_category;
	
	I3Matrix() {};
	I3Matrix(const base &m) : base(m) {};
	I3Matrix(size_t size1, size_t size2) : base(size1, size2) {}
	I3Matrix(size_t size1, size_t size2, const value_type &init) : base(size1, size2, init) {}
	I3Matrix(size_t size1, size_t size2, const array_type &data) : base(size1, size2, data) {}

private:
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

I3_POINTER_TYPEDEFS(I3Matrix);
// BOOST_CLASS_VERSION(I3Matrix, 0);

#endif // I3CLSIM_I3MATRIX_H_INCLUDED
