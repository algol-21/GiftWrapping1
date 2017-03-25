#include "MathVector.h"

MathVector::MathVector()
{
	math_vector = std::vector<double>();
}

MathVector::MathVector(const std::vector<double>& _math_vector_)
{
	math_vector = _math_vector_;
}

MathVector::MathVector(size_t dimension)
{
	math_vector = std::vector<double>(dimension);
}

MathVector::~MathVector()
{
}

size_t MathVector::getDimension() const 
{
	return math_vector.size();
}

void MathVector::normalize()
{
	double length = sqrt((*this) * (*this));
	*this = (1.0 / length) * (*this);
}

MathVector operator+(const MathVector& left, const MathVector& right)
{
	MathVector result = left;
	for (size_t counter = 0; counter < result.math_vector.size(); ++counter)
		result.math_vector[counter] += right.math_vector[counter];

	return result;
}

MathVector operator*(double scalar, const MathVector& vector)
{
	MathVector result = vector;
	for (size_t counter = 0; counter < result.math_vector.size(); ++counter)
		result.math_vector[counter] *= scalar;

	return result;
}

MathVector operator-(const MathVector& left, const MathVector& right)
{
	 return left + (-1)*right;
}

double operator*(const MathVector& left, const MathVector& right)
{
	double result = 0.0;
	for (size_t counter = 0; counter < left.math_vector.size(); ++counter)
		result += left.math_vector[counter] * right.math_vector[counter];

	return result;
}

std::ostream& operator<<(std::ostream& output_stream, const MathVector& vector)
{
	//output_stream << "{ ";
	for (auto coordinate : vector.math_vector)
		output_stream << coordinate << " ";
	//output_stream << "}";

	return output_stream;
}

double& MathVector::operator[](size_t index)
{
	return math_vector[index];
}

double MathVector::operator[](size_t index) const
{
	return math_vector[index];
}

MathVector MathVector::crossProduct(const std::vector<MathVector>& vectors)
{
	size_t dimension = vectors[0].getDimension();
	MathVector result(dimension);
	std::vector<std::vector<double>> tmp_matrix(dimension - 1, std::vector<double>(dimension - 1));

	size_t row_tmp_matrix;
	bool flag_skip;

	for (size_t index = 0; index < dimension; ++index)
	{
		flag_skip = false;

		for (size_t row = 0; row < dimension; ++row)
		{
			if (row == index)
			{
				flag_skip = true;
				continue;
			}

			row_tmp_matrix = row;
			if (flag_skip)
				--row_tmp_matrix;
				

			for (size_t column = 0; column < dimension - 1; ++column)
				tmp_matrix[column][row_tmp_matrix] = vectors[column][row];
		}

		result[index] = pow(-1, index + 1 + dimension) * determinant(tmp_matrix);
	}

	return result;
}

// Передаём ссылку для того, чтобы дополнительно не вызывать конструктор копирования.
// Это не совсем корректно, т.к. изменяется внешний объект matrix, но determinant используется только в crossProduct,
// а после каждого вызова determinant matrix перед последующим использованием полностью пересобирается.
double MathVector::determinant(std::vector<std::vector<double>>& matrix)
{
	double result = 1.0;
	double tmp;
	bool flag_NonZeroExits;

	for (size_t row_1 = 0; row_1 < matrix.size(); ++row_1)
	{
		if (matrix[row_1][row_1] == 0)
		{
			flag_NonZeroExits = false;

			for (size_t row_2 = row_1 + 1; row_2 < matrix.size(); ++row_2)
				if (matrix[row_2][row_1] != 0)
				{
					for (size_t column = row_1; column < matrix.size(); ++column)
					{
						tmp = matrix[row_1][column];
						matrix[row_1][column] = matrix[row_2][column];
						matrix[row_2][column] = tmp;
						//result *= -1.0;
						//flag_NonZeroExits = true;
					}

					// -----
					result *= -1.0;
					flag_NonZeroExits = true;
					// -----
					break;
				}
			
			if (!flag_NonZeroExits)
				return 0.0;
		}

		for (size_t row_2 = row_1 + 1; row_2 < matrix.size(); ++row_2)
			for (size_t column = row_1 + 1; column < matrix.size(); ++column)
				matrix[row_2][column] -= matrix[row_1][column] * matrix[row_2][row_1] / matrix[row_1][row_1];
	}

	for (size_t index = 0; index < matrix.size(); ++index)
		result *= matrix[index][index];

	return result;
}