#include <iostream>
#include <vector>
#include "iomanip"
#include "cmath"

using namespace std;

const int INF = 2e9;

class Matrix {
protected:
    int n, m;
    vector<vector<double>> matrix;
public:
    Matrix() {
        n = 0;
        m = 0;
    }

    Matrix(int rows, int cols, vector<vector<double>> &mas) {
        n = rows;
        m = cols;
        for (int i = 0; i < n; ++i) {
            matrix.push_back(vector<double>(m));
            for (int j = 0; j < m; ++j) {
                matrix[i][j] = mas[i][j];
            }
        }
    }

    Matrix(int rows, int cols, vector<double> &arr) {
        n = rows;
        m = cols;
        for (int i = 0; i < n; ++i) {
            matrix.push_back(vector<double>(m));
            for (int j = 0; j < m; ++j) {
                matrix[i][j] = pow(arr[i], j);
            }
        }
    }

    Matrix(Matrix &a, Matrix &b) {
        matrix.clear();
        n = a.getRows();
        m = a.getColumns() + b.getColumns();
        for (int i = 0; i < n; ++i) {
            matrix.push_back(vector<double>(m));
            for (int j = 0; j < a.getColumns(); ++j) {
                matrix[i][j] = a.getElement(i, j);
            }
            for (int j = 0; j < b.getColumns(); ++j) {
                matrix[i][a.getColumns() + j] = b.getElement(i, j);
            }
        }
    }

    virtual int getMaxAbsRow(int column) {
        double maximum = getElement(column, column);
        int index = column;
        for (int i = column; i < n; ++i) {
            if (abs(maximum) < abs(getElement(i, column))) {
                index = i;
                maximum = getElement(i, column);
            }
        }
        return index;
    }

    virtual int getRows() {
        return n;
    }

    virtual int getColumns() {
        return m;
    }

    virtual double getElement(int i, int j) {
        if (i >= n || j >= m || j < 0 || i < 0) {
            cout << "Error: the dimensional problem occurred\n";
            return -INF;
        }
        return matrix[i][j];
    }

    virtual Matrix operator+(Matrix other) {
        vector<vector<double>> temp;
        if (n != other.getRows() || m != other.getColumns()) {
            cout << "Error: the dimensional problem occurred\n";
            return Matrix();
        }
        for (int i = 0; i < n; ++i) {
            temp.push_back(vector<double>(m));
            for (int j = 0; j < m; ++j)
                temp[i][j] = matrix[i][j] + other.getElement(i, j);
        }
        return Matrix(n, m, temp);
    }

    virtual Matrix operator-(Matrix other) {
        vector<vector<double>> temp;
        if (n != other.getRows() || m != other.getColumns()) {
            cout << "Error: the dimensional problem occurred\n";
            return Matrix();
        }
        for (int i = 0; i < n; ++i) {
            temp.push_back(vector<double>(m));
            for (int j = 0; j < m; ++j)
                temp[i][j] = matrix[i][j] - other.getElement(i, j);
        }
        return Matrix(n, m, temp);
    }

    virtual Matrix operator*(Matrix other) {
        double sum;
        vector<vector<double>> temp;
        if (m != other.getRows()) {
            cout << "Error: the dimensional problem occurred\n";
            return Matrix();
        }
        for (int i = 0; i < n; ++i) {
            temp.push_back(vector<double>(m));
            for (int j = 0; j < other.getColumns(); ++j) {
                sum = 0;
                for (int p = 0; p < m; ++p)
                    sum += matrix[i][p] * other.getElement(p, j);
                temp[i][j] = sum;
            }
        }
        return Matrix(n, other.getColumns(), temp);
    }

    void operator=(Matrix other) {
        matrix.resize(other.getRows());
        for (int i = 0; i < other.getRows(); ++i) {
            matrix[i].resize(other.getColumns());
            for (int j = 0; j < other.getColumns(); ++j)
                matrix[i][j] = other.getElement(i, j);
        }
        n = other.getRows();
        m = other.getColumns();
    }

    virtual Matrix transpose() {
        vector<vector<double>> temp;
        for (int i = 0; i < m; ++i) {
            temp.push_back(vector<double>(n));
            for (int j = 0; j < n; ++j) {
                temp[i][j] = matrix[j][i];
            }
        }
        return Matrix(m, n, temp);
    }
};

class ColumnVector : public Matrix {
public:
    ColumnVector() = default;

    ColumnVector(int s, vector<double> elements) {
        vector<vector<double>> temp;
        for (int i = 0; i < s; ++i) {
            temp.push_back(vector<double>(1, elements[i]));
        }
        m = 1;
        n = s;
        Matrix *self = this;
        *self = Matrix(n, m, temp);
    }

    double norm() {
        double sum = 0;
        for (int i = 0; i < n; ++i)
            sum += pow(matrix[i][0], 2.0);
        return sqrt(sum);
    }

    Matrix operator+(ColumnVector other) {
        Matrix *my = this;
        Matrix *arg = &other;
        return *my + *arg;
    }

    Matrix operator-(ColumnVector other) {
        Matrix *my = this;
        Matrix *arg = &other;
        return *my - *arg;
    }

    Matrix transpose() override {
        Matrix *my = this;
        return my->Matrix::transpose();
    }

    bool operator==(Matrix other) {
        if (n != other.getRows() || m != other.getColumns())
            return false;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                if (matrix[i][j] != other.getElement(i, j))
                    return false;
        return true;
    }

    void operator=(Matrix other) {
        Matrix *my = this;
        *my = other;
    }

    Matrix operator*(ColumnVector other) {
        Matrix *my = this;
        Matrix *arg = &other;
        return *my * *arg;
    }
};

class SquareMatrix : public Matrix {
public:
    SquareMatrix() = default;

    SquareMatrix(int e, vector<vector<double>> &mas) {
        n = e;
        m = e;
        for (int i = 0; i < n; ++i) {
            matrix.push_back(vector<double>(m));
            for (int j = 0; j < m; ++j) {
                matrix[i][j] = mas[i][j];
            }
        }
    }

    SquareMatrix(Matrix a) {
        Matrix *my = this;
        *my = a;
    }

    Matrix operator+(SquareMatrix other) {
        Matrix *my = this;
        Matrix *arg = &other;
        return *my + *arg;
    }

    Matrix operator-(SquareMatrix other) {
        Matrix *my = this;
        Matrix *arg = &other;
        return *my - *arg;
    }

    Matrix transpose() override {
        Matrix *my = this;
        return my->Matrix::transpose();
    }

    bool operator==(Matrix other) {
        if (n != other.getRows() || m != other.getColumns())
            return false;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                if (matrix[i][j] != other.getElement(i, j))
                    return false;
        return true;
    }

    void operator=(Matrix other) {
        Matrix *my = this;
        *my = other;
    }

    Matrix operator*(SquareMatrix other) {
        Matrix *my = this;
        Matrix *arg = &other;
        return *my * *arg;
    }
};

class IdentityMatrix : public SquareMatrix {
public:
    IdentityMatrix() = default;

    IdentityMatrix(int elements) {
        vector<vector<double>> temp;
        for (int i = 0; i < elements; ++i) {
            temp.push_back(vector<double>(elements, 0));
            temp[i][i] = 1;
        }
        SquareMatrix *self = this;
        *self = SquareMatrix(elements, temp);
    }
};

class DiagonalInverseMatrix : public SquareMatrix {
public:
    DiagonalInverseMatrix(int elements, Matrix a) {
        vector<vector<double>> temp;
        for (int i = 0; i < elements; ++i) {
            temp.push_back(vector<double>(elements, 0));
            temp[i][i] = 1 / a.getElement(i, i);
        }
        SquareMatrix *self = this;
        *self = SquareMatrix(elements, temp);
    }
};

class EliminationMatrix : public IdentityMatrix {
public:
    EliminationMatrix() = default;

    EliminationMatrix(int elements, SquareMatrix a, int i, int j) {
        IdentityMatrix identity = IdentityMatrix(elements);
        Matrix *casted = &identity;
        Matrix *self = this;
        *self = *casted;
        this->matrix[i][j] = -a.getElement(i, j) / a.getElement(j, j);
    }
};

class PermutationMatrix : public IdentityMatrix {
public:
    PermutationMatrix() = default;

    PermutationMatrix(int elements, int i, int j) {
        IdentityMatrix identity = IdentityMatrix(elements);
        Matrix *self = this;
        Matrix *casted = &identity;
        *self = *casted;
        this->matrix[i][i] = 0;
        this->matrix[j][j] = 0;
        this->matrix[i][j] = 1;
        this->matrix[j][i] = 1;
    }
};

ostream &operator<<(ostream &out, Matrix matr) {
    for (int i = 0; i < matr.getRows(); ++i) {
        for (int j = 0; j < matr.getColumns(); ++j) {
            if (abs(matr.getElement(i, j)) < 0.0001)
                out << "0.0000 ";
            else
                out << fixed << setprecision(4) << matr.getElement(i, j) << ' ';
        }
        out << '\n';
    }
    return out;
}

istream &operator>>(istream &is, Matrix &matr) {
    int n, m;
    vector<vector<double>> temp;
    is >> n >> m;
    for (int i = 0; i < n; ++i) {
        temp.push_back(vector<double>(m));
        for (int j = 0; j < m; ++j) {
            is >> temp[i][j];
        }
    }
    matr = Matrix(n, m, temp);
    return is;
}

istream &operator>>(istream &in, SquareMatrix &matr) {
    int n;
    vector<vector<double>> temp;
    in >> n;
    for (int i = 0; i < n; ++i) {
        temp.push_back(vector<double>(n));
        for (int j = 0; j < n; ++j) {
            in >> temp[i][j];
        }
    }
    matr = SquareMatrix(n, temp);
    return in;
}

ostream &operator<<(ostream &os, SquareMatrix matr) {
    Matrix *arg = &matr;
    os << *arg;
    return os;
}

ostream &operator<<(ostream &os, IdentityMatrix matr) {
    Matrix *arg = &matr;
    os << *arg;
    return os;
}

ostream &operator<<(ostream &os, EliminationMatrix matr) {
    Matrix *arg = &matr;
    os << *arg;
    return os;
}

ostream &operator<<(ostream &os, PermutationMatrix matr) {
    Matrix *arg = &matr;
    os << *arg;
    return os;
}

void inverse(SquareMatrix &a) {
    int step = 1;
    IdentityMatrix identity = IdentityMatrix(a.getColumns());
    for (int i = 0; i < a.getColumns(); ++i) {
        int index = a.getMaxAbsRow(i);
        if (index > i) {
            PermutationMatrix p = PermutationMatrix(a.getColumns(), index, i);
            a = p * a;
            *(Matrix *) &identity = p * identity;
            ++step;
        }
        for (int j = i + 1; j < a.getRows(); ++j) {
            EliminationMatrix e = EliminationMatrix(a.getColumns(), a, j, i);
            if (a.getElement(j, i) == 0)
                continue;
            a = e * a;
            *(Matrix *) &identity = e * identity;
            ++step;
        }
    }
    for (int i = a.getColumns() - 1; i >= 0; --i) {
        for (int j = i - 1; j >= 0; --j) {
            if (a.getElement(j, i) == 0)
                continue;
            EliminationMatrix e = EliminationMatrix(a.getColumns(), a, j, i);
            a = e * a;
            *(Matrix *) &identity = e * identity;
            ++step;
        }
    }

    DiagonalInverseMatrix inverse = DiagonalInverseMatrix(a.getColumns(), a);
    a = inverse * a;
    *(Matrix *) &identity = inverse * identity;
    a = identity;
}

int main() {
    int n, degree;
    vector<double> datX;
    vector<double> datY;
    double x, y;
    cin >> n;
    for (int i = 0; i < n; ++i) {
        cin >> x >> y;
        datX.push_back(x);
        datY.push_back(y);
    }
    cin >> degree;
    Matrix a = Matrix(n, degree + 1, datX);
    ColumnVector b = ColumnVector(n, datY);
    cout << "A:\n" << a;
    SquareMatrix multipliedByTranspose = a.transpose() * a;
    cout << "A_T*A:\n" << multipliedByTranspose;
    inverse(multipliedByTranspose);
    cout << "(A_T*A)^-1:\n";
    cout << multipliedByTranspose;
    ColumnVector transposedByB;
    transposedByB = a.transpose() * b;
    cout << "A_T*b:\n" << transposedByB;
    cout << "x~:\n" << multipliedByTranspose * transposedByB;
    return 0;
}
