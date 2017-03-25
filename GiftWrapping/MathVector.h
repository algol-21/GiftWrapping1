#pragma once

#include <vector>
#include <iostream>

class MathVector
{
public:
	MathVector();
	MathVector(const std::vector<double>&);
	MathVector(size_t);
	~MathVector();

	size_t getDimension() const;
	void normalize();
	static MathVector crossProduct(const std::vector<MathVector>&);

	friend MathVector    operator+(const MathVector&, const MathVector&);
	friend MathVector    operator-(const MathVector&, const MathVector&);
	friend MathVector    operator*(double, const MathVector&);
	friend double        operator*(const MathVector&, const MathVector&);
	friend std::ostream& operator<<(std::ostream&, const MathVector&);
	
	double& operator[](size_t);
	double  operator[](size_t) const;

private:
	std::vector<double> math_vector;

	static double determinant(std::vector<std::vector<double>>&);
};

