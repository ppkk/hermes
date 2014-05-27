#ifndef IML_WRAPPERS_H
#define IML_WRAPPERS_H

#include "hermes_common.h"
#include "hermes2d.h"
#include <iostream>

#include "cg.h"
#include "gmres.h"

class IMLVector
{
public:
    IMLVector()
    {
        m_size = 0;
        m_data = nullptr;
    }

    explicit IMLVector(int n)
    {
        m_size = 0;
        m_data = nullptr;
        alloc(n);
    }

    IMLVector(const IMLVector& vec)
    {
        m_size = 0;
        m_data = nullptr;
        alloc(vec.size());
        memcpy(this->m_data, vec.m_data, size() * sizeof(double));
    }

    ~IMLVector()
    {
        free();
    }

    void free()
    {
        if(m_data)
            delete[] m_data;
        m_data = nullptr;
        m_size = 0;
    }

    void alloc(int n)
    {
        free();
        m_size = n;
        m_data = new double[n];
        memset(m_data, 0, n*sizeof(double));
    }

    void set(SimpleVector<double>* x)
    {
        free();
        alloc(x->get_size());
        memcpy(m_data, x->v, m_size * sizeof(double));
    }

    void set(double* source, int size)
    {
        free();
        alloc(size);
        memcpy(m_data, source, m_size * sizeof(double));
    }

    int size() const
    {
        return m_size;
    }

    void copy(SimpleVector<double>* dest) const
    {
        dest->alloc(m_size);
        dest->set_vector(m_data);
    }

    void operator= (const IMLVector& x)
    {
        if(m_size == 0)
            alloc(x.size());
        assert(m_size == x.size());

        (memcpy(m_data, x.m_data, m_size * sizeof(double)));
    }

    void operator= (double x)
    {
        for(int i = 0; i < size(); i++)
            (*this)(i) = x;
    }

    double& operator() (int idx)
    {
        return m_data[idx];
    }

    double operator() (int idx) const
    {
        return m_data[idx];
    }

    IMLVector operator-(const IMLVector& b) const
    {
        if(b.size() != this->size())
        {
            std::cout << "operator- : incompatible sizes " << b.size() << " !=  " << this->size() << std::endl;
            assert(0);
        }

        IMLVector result;
        result.alloc(size());
        for(int i = 0; i < size(); i++)
            result(i) = (*this)(i) - b(i);

        return result;
    }

    IMLVector operator+(const IMLVector& b) const
    {
        if(b.size() != this->size())
        {
            std::cout << "operator- : incompatible sizes " << b.size() << " !=  " << this->size() << std::endl;
            assert(0);
        }

        IMLVector result;
        result.alloc(size());
        for(int i = 0; i < size(); i++)
            result(i) = (*this)(i) + b(i);

        return result;
    }

    void operator-=(const IMLVector& b)
    {
        assert(b.size() == this->size());
        for(int i = 0; i < size(); i++)
            (*this)(i) -= b(i);

    }

    void operator+=(const IMLVector& b)
    {
        assert(b.size() == this->size());
        for(int i = 0; i < size(); i++)
            (*this)(i) += b(i);

    }

    IMLVector operator*(double alpha) const
    {
        IMLVector result;
        result.alloc(size());
        for(int i = 0; i < size(); i++)
            result(i) = (*this)(i) * alpha;

        return result;
    }

    friend std::ostream& operator<<(std::ostream& os, const IMLVector& vec);

    double* m_data;
    int m_size;
};

std::ostream& operator<<(std::ostream& os, const IMLVector& vec)
{
    os << "[";
    for(int i = 0; i < vec.size(); i++)
    {
        os << vec(i);
        if(i < vec.size() - 1)
            os << ", ";
    }
    os << "]";
    return os;
}


IMLVector operator*(double alpha, const IMLVector& vec)
{
    IMLVector result;
    result.alloc(vec.size());
    for(int i = 0; i < vec.size(); i++)
        result(i) = vec(i) * alpha;

    return result;
}

double dot(const IMLVector& a, const IMLVector& b)
{
    assert(a.size() == b.size());
    double res = 0;

    for(int i = 0; i < a.size(); i++)
        res += a(i) * b(i);

    return res;
}

double norm(const IMLVector& x)
{
    return std::sqrt(dot(x, x));
}

class IMLOperator
{
public:
    IMLOperator(CSCMatrix<double>* csc_matrix) : m_csc_matrix(csc_matrix) {}

    IMLVector operator*(const IMLVector& x) const
    {
        assert(x.size() == m_csc_matrix->get_size());
        SimpleVector<double> result, x_simple_vec;
        result.alloc(x.size());
        x_simple_vec.alloc(x.size());
        x.copy(&x_simple_vec);
        m_csc_matrix->multiply_with_vector(&x_simple_vec, &result);
        IMLVector result_iml;
        result_iml.set(&result);
        return result_iml;
    }

    CSCMatrix<double>* m_csc_matrix;
};



class IMLEmptyPreconditioner
{
public:

    IMLVector solve(const IMLVector& x) const
    {
        return x;
    }
};

class IMLDiagPreconditioner
{
public:
    IMLDiagPreconditioner(CSCMatrix<double>* matrix) : matrix(matrix) {}

    IMLVector solve(const IMLVector& x) const
    {
        assert(x.size() == matrix->get_size());
        SimpleVector<double> diagonal;
        diagonal.alloc(matrix->get_size());
        matrix->extract_diagonal(&diagonal);
        //m_matrix->export_to_file("matrix.m", "Matrix", EXPORT_FORMAT_PLAIN_ASCII);
        IMLVector result;
        result.alloc(matrix->get_size());
        for(int i = 0; i < matrix->get_size(); i++)
        {
            result(i) = x(i) / diagonal.v[i];
        }

        return result;
    }

    CSCMatrix<double>* matrix;
};


class IMLMatrix
{
public:
    IMLMatrix(int n, int m) : nrows(n), ncols(m)
    {
        data = new double[n*m];
        memset(data, 0, n*m*sizeof(double));

        matrix_changed = true;
        umfpack_matrix = nullptr;
        umfpack_solver = nullptr;
        umfpack_rhs = nullptr;
    }
    ~IMLMatrix()
    {
        delete[] data;

        if(umfpack_matrix)
            delete umfpack_matrix;
        if(umfpack_rhs)
            delete umfpack_rhs;
        if(umfpack_solver)
            delete umfpack_solver;
    }

    IMLMatrix transpose()
    {
        IMLMatrix result(ncols, nrows);
        for(int col = 0; col < ncols; col++)
        {
            for(int row = 0; row < nrows; row++)
            {
                result(col, row) = (*this)(row, col);
            }
        }

        return result;
    }

    double& operator() (int i, int j)
    {
        assert( (i >= 0) && (i < nrows) && (j >= 0) && (j < ncols) );
        matrix_changed = true;
        return data[ncols*i + j];
    }

    double operator() (int i, int j) const
    {
        assert( (i >= 0) && (i < nrows) && (j >= 0) && (j < ncols) );
        return data[ncols*i + j];
    }

    IMLVector operator*(const IMLVector& x) const
    {
        assert(x.size() == ncols);
        IMLVector result(nrows);
        for(int i = 0; i < nrows; i++)
        {
            result(i) = 0.;
            for(int j = 0; j < ncols; j++)
            {
                result(i) += (*this)(i,j) * x(j);
            }
        }
        return result;
    }

    IMLVector solve(const IMLVector& arg)
    {
        if(!umfpack_solver)
        {
            assert(!umfpack_matrix);
            assert(!umfpack_rhs);
            umfpack_matrix = new CSCMatrix<double>;
            umfpack_rhs = new SimpleVector<double>;
            umfpack_solver = new Hermes::Solvers::UMFPackLinearMatrixSolver<double>(umfpack_matrix, umfpack_rhs);
        }

        arg.copy(umfpack_rhs);
        if(matrix_changed)
        {
            umfpack_matrix->zero();
            assert(nrows == ncols);
            int size = nrows;
            int nnz = size*size;
            int* ap = new int[size+1];
            int* ai = new int[nnz];
            double* ax = new double[nnz];

            for(int i = 0; i <= size; i++)
            {
                ap[i] = i * size;
            }
            for(int i = 0; i < size; i++)
            {
                for(int j = 0; j < size; j++)
                {
                    ai[size*i + j] = j;
                    ax[size*i + j] = (*this)(j,i);
                }
            }

            umfpack_matrix->create(size, nnz, ap, ai, ax);
        }


        if(matrix_changed)
            umfpack_solver->set_reuse_scheme(Hermes::Solvers::HERMES_CREATE_STRUCTURE_FROM_SCRATCH);
        else
            umfpack_solver->set_reuse_scheme(Hermes::Solvers::HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);

        umfpack_solver->solve();
        matrix_changed = false;

        IMLVector result;
        result.set(umfpack_solver->get_sln_vector(), nrows);

        return result;
    }

    IMLVector solve_gmres(const IMLVector& arg)
    {
        IMLEmptyPreconditioner empty_preconditioner;
        int gmres_m = 50;
        double tol = 1e-12;
        int max_iter = 500;
        int size = arg.size();
        IMLMatrix iml_matrix(gmres_m+1, size);
        IMLVector iml_x(size);

        int converged = GMRES(*this, iml_x, arg, empty_preconditioner, iml_matrix, gmres_m, max_iter, tol);
        Hermes::Mixins::Loggable::Static::info("GMRES converged %d, tol %g, steps %d, ndofs %d", converged, tol, max_iter, size);

        return iml_x;
    }

    void print_sparse(const char* filename, const char* varname)
    {
      return;
      int* rows = new int[nrows*ncols];
      int* cols = new int[nrows*ncols];
      double* vals = new double[nrows*ncols];

        int nnz = 0;
        for(int row = 0; row < nrows; row++)
        {
            for(int col = 0; col < ncols; col++)
            {
                if(fabs((*this)(row, col)) > 1e-13)
                {
                    rows[nnz] = row+1;
                    cols[nnz] = col+1;
                    vals[nnz] = (*this)(row, col);
                    nnz++;
                }
            }
        }

        FILE* f = NULL;
        f = fopen(filename, "w");
        assert(f!=NULL);

        fprintf(f, "rows=[");
        for(int i = 0; i < nnz; i++)
            fprintf(f, "%d, ", rows[i]);
        fprintf(f, "];\ncols=[");
        for(int i = 0; i < nnz; i++)
            fprintf(f, "%d, ", cols[i]);
        fprintf(f, "];\nvals = [");
        for(int i = 0; i < nnz; i++)
            fprintf(f, "%g, ", vals[i]);
        fprintf(f, "];\n %s = sparse(rows, cols, vals);", varname);
        fclose(f);
    }

    friend std::ostream& operator<<(std::ostream& os, const IMLMatrix& mat);

    int num_rows() const { return nrows; }
    int num_cols() const { return ncols; }

    double *data;
    int nrows, ncols;

    // used when solve is required
    Solvers::UMFPackLinearMatrixSolver<double> *umfpack_solver;
    CSCMatrix<double> *umfpack_matrix;
    SimpleVector<double> *umfpack_rhs;
    bool matrix_changed;
};

std::ostream& operator<<(std::ostream& os, const IMLMatrix& mat)
{
    os << "*********************************" << std::endl;
    for(int i = 0; i < mat.nrows; i++)
    {
        for(int j = 0; j < mat.ncols; j++)
        {
            os << mat(i, j) << " ";
        }
        os <<";" << std::endl;
    }
    os << "*********************************" << std::endl;
    return os;
}


template <class Operator, class Precond>
class IMLPrecondTimesMatrixOperator
{
public:
    IMLPrecondTimesMatrixOperator(const Operator* oper, const Precond* precond) : oper(oper), precond(precond) {}
    IMLVector operator*(const IMLVector& x) const
    {
        IMLVector result;
        result = (*precond).solve((*oper)*x);
        return result;
    }

    const Operator* oper;
    const Precond* precond;
};

template <class Operator, class Matrix>
void IMLOperatorToMatrix(const Operator& oper, Matrix &mat)
{
    for(int col = 0; col < mat.num_cols(); col++)
    {
        IMLVector vec(mat.num_cols());
        vec(col) = 1.;
        IMLVector oper_times_vec(mat.num_rows());
        oper_times_vec = oper * vec;
        for(int row = 0; row < mat.num_rows(); row++)
        {
            mat(row, col) = oper_times_vec(row);
        }
    }
}

void test_iml_wrappers()
{
    int ndofs = 4;
    IMLMatrix mat(ndofs, ndofs);

    mat(0,0) = 5;
    mat(1,1) = 9;
    mat(2,2) = 44;
    mat(3,3) = 1;
    mat(3,2) = -5;
    mat(1,3) = 7;
    mat(2, 3) = 1;
    mat(0,2) = -1;
    std::cout << mat << std::endl;
    IMLVector vec(ndofs), x(ndofs);
    vec(0) = 23;
    vec(1) = -44;
    vec(2) = 11;
    //std::cout << mat*vec << std::endl;

    IMLEmptyPreconditioner empty_preconditioner;
    double tol = 1e-6;
    int max_iter = 50;
    int gmres_m = 10;
    IMLMatrix helper_matrix(50, 50);
    int converged = GMRES(mat,  x, vec, empty_preconditioner, helper_matrix, gmres_m, max_iter, tol);
    //int converged = CG(mat,  x, vec, iml_preconditioner, max_iter, tol);
    std::cout << "converged " << converged << ", tol " << tol << ", steps " << max_iter << std::endl;
    std::cout << "rhs   " << vec << std::endl;
    std::cout << "x     " << x << std::endl;
    std::cout << "mat*x " << mat * x << std::endl;

    std::cout << "Testing IMLOperatorToMatrix:\n";
    IMLMatrix mat_cpy(ndofs, ndofs);
    std::cout << mat;
    IMLOperatorToMatrix(mat, mat_cpy);
    std::cout << mat_cpy;

    x = 0.0;
    std::cout << "x     " << x << std::endl;
    std::cout << "rhs   " << vec << std::endl;
    x = mat_cpy.solve(vec);
    std::cout << "x     " << x << std::endl;
    std::cout << "mat*x " << mat * x << std::endl;

    x = 0.0;
    std::cout << "x     " << x << std::endl;
    std::cout << "rhs   " << vec << std::endl;
    x = mat_cpy.solve_gmres(vec);
    std::cout << "x     " << x << std::endl;
    std::cout << "mat*x " << mat * x << std::endl;


}

#endif // IML_WRAPPERS_H
