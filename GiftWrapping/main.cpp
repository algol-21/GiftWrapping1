#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <unordered_set>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/math/constants/constants.hpp>
#include "MathVector.h"

size_t us_hash(const std::unordered_set<size_t>& S)
{
	size_t sum = 0;
	for (auto item : S)
		sum += item;
	return sum;
}

// void find_subfaces(const std::vector<size_t>& hyperface, std::vector<std::vector<size_t>>& subfaces)
void find_subfaces(const std::unordered_set<size_t>& hyperface, std::unordered_set<std::unordered_set<size_t>, decltype(&us_hash)>& subfaces)
{
// -----
	subfaces.clear();
// -----

	// -----
	// std::vector<size_t> subface;
	std::unordered_set<size_t> subface;
	// -----


	for (size_t indexSkip = 0; indexSkip < hyperface.size(); ++indexSkip)
	{
		subface.clear();

		// -----
		auto it = hyperface.begin();
		// -----

		for (size_t index = 0; index < hyperface.size(); ++index)
		{
			if (index != indexSkip)
				// -----
				// subface.push_back(hyperface[index]);
				subface.insert(*it);
				// -----
			++it;
		}

		// -----
		// subfaces.push_back(subface);
		subfaces.insert(subface);
		// -----
	}
}

// -----
// all_points изменяется в несимлициальном случае, потому std::vector<MathVector>& all_points (без const)
size_t wrapping(bool is_first_hyperface, const MathVector& normal_of_hyperface, const MathVector& point_of_subface, const MathVector& normal_of_subface, const std::unordered_set<size_t>& indexes_of_candidates, std::vector<MathVector>& all_points)
{
	// ~~~~~
	double eps = 1e-10;
	// ~~~~~

	MathVector vector_in_new_candidate_hyperface = all_points[*indexes_of_candidates.begin()] - point_of_subface;
	
	// ~~~~~
	// =====
	double dot_v_a = vector_in_new_candidate_hyperface * normal_of_subface;
	// =====
	double abs_dot_v_n = fabs(vector_in_new_candidate_hyperface * normal_of_hyperface);

	if (abs_dot_v_n < eps)
	{
		if (!is_first_hyperface)
			std::cout << "A non-simplicial case is found" << std::endl;
		// =====
		//if (dot_v_a > 0)
			return *indexes_of_candidates.begin();
		//else
		//	std::cout << "wrapping warning: \"ctg\" should be = +inf" << std::endl;
		// =====
	}
		
	// ~~~~~

	double ctg = -dot_v_a / abs_dot_v_n;
	double min_ctg = ctg;
	size_t index_of_point_with_min_ctg = *indexes_of_candidates.begin();

	for (auto it = std::next(indexes_of_candidates.begin(), 1); it != indexes_of_candidates.end(); ++it)
	{
		vector_in_new_candidate_hyperface = all_points[*it] - point_of_subface;
		
		// ~~~~~

		// =====
		dot_v_a = vector_in_new_candidate_hyperface * normal_of_subface;
		// =====

		abs_dot_v_n = fabs(vector_in_new_candidate_hyperface * normal_of_hyperface);

		if (abs_dot_v_n < eps)
		{
			if (!is_first_hyperface)
			{
				std::cout << "A non-simplicial case is found" << std::endl;

				// To generate a random offset vector ...
				
				double radius = 1;//1e-3;
				size_t dimension = all_points[0].getDimension();

				boost::mt19937 generator;
				auto distribution = boost::uniform_real<>(0.0, 2.0 * boost::math::constants::pi<double>());

				std::vector<double> angles(dimension - 1);
				
				for (size_t i = 0; i < dimension - 1; ++i)
					angles[i] = distribution(generator);


				//std::vector<double> offset(dimension, radius);
				MathVector offset(std::vector<double>(dimension, radius));

				for (size_t coordinate = 0; coordinate < dimension; ++coordinate)
				{
					for (size_t i = 0; i < dimension - coordinate - 1; ++i)
						offset[coordinate] *= cos(angles[i]);

					if (coordinate != 0)
						offset[coordinate] *= sin(angles[dimension - coordinate - 1]);
				}

				// Get a moved point
				// -----
				//all_points[*it] += offset;
				all_points[*it] = all_points[*it] + offset;
				// -----
			}
				
			// =====
			//if (dot_v_a > 0)
				return *it;
			//else
			//	continue;
			// =====
		}

		ctg = -dot_v_a / abs_dot_v_n;
		// ~~~~~

		if (ctg < min_ctg)
		{
			min_ctg = ctg;
			index_of_point_with_min_ctg = *it;
		}
 
	}

	return index_of_point_with_min_ctg;
}


void create_coordinate_axis(size_t num_of_coordinate, size_t dimension, MathVector& coordinate_axis)
{
	coordinate_axis = std::vector<double>(dimension, 0.0);
	coordinate_axis[num_of_coordinate] = 1.0;
}

// Формально меняет all_points, но на самом деле нет
void find_first_hyperface(std::vector<MathVector>& all_points, std::unordered_set<size_t>& first_hyperface)
{
	// Initialize remaining indexes, which can be considered.
	std::unordered_set<size_t> indexes_of_candidates;
	for (size_t counter = 0; counter < all_points.size(); ++counter)
		indexes_of_candidates.insert(counter);

	// Find point with min first coordinate. First normal is (1, 0, ..., 0).
	double min_first_coordinate = all_points[0][0];
	size_t index_of_point_with_min_first_coordinate = 0;

	for (size_t counter = 1; counter < all_points.size(); ++counter)
		if (all_points[counter][0] < min_first_coordinate)
		{
			min_first_coordinate = all_points[counter][0];
			index_of_point_with_min_first_coordinate = counter;
		}

	first_hyperface.insert(index_of_point_with_min_first_coordinate);
	indexes_of_candidates.erase(index_of_point_with_min_first_coordinate);


	MathVector normal_of_subface;
	MathVector normal_of_hyperface;

	size_t dimension = all_points[0].getDimension();
	create_coordinate_axis(0, dimension, normal_of_hyperface);
	
	
	std::vector<MathVector> cross_product_vectors;
	MathVector coordinate_axis;
	size_t new_index;

	for (size_t counter_1 = 1; counter_1 < dimension; ++counter_1)
	{
		// Recalculate normal of a subface
		cross_product_vectors.clear();
		cross_product_vectors.push_back(normal_of_hyperface);

		auto it = std::next(first_hyperface.begin());

		for (size_t counter_2 = 1; counter_2 < counter_1; ++counter_2)
		{
			cross_product_vectors.push_back(all_points[*it] - all_points[*first_hyperface.begin()]);
			++it;
		}

		for (size_t counter_2 = counter_1 + 1; counter_2 < dimension; ++counter_2)
		{
			create_coordinate_axis(counter_2, dimension, coordinate_axis);
			cross_product_vectors.push_back(coordinate_axis);
		}

		normal_of_subface = MathVector::crossProduct(cross_product_vectors);
		normal_of_subface.normalize();

		new_index = wrapping(true, normal_of_hyperface, all_points[*first_hyperface.begin()], normal_of_subface, indexes_of_candidates, all_points);
		first_hyperface.insert(new_index);
		indexes_of_candidates.erase(new_index);

		// Recalculate normal of a hyperface
		cross_product_vectors.clear();
	
		it = std::next(first_hyperface.begin(), 1);
		for (size_t counter_2 = 1; counter_2 < counter_1 + 1; ++counter_2)
		{
			cross_product_vectors.push_back(all_points[*it] - all_points[*first_hyperface.begin()]);
			++it;
		}

		for (size_t counter_2 = counter_1 + 1; counter_2 < dimension; ++counter_2)
		{
			create_coordinate_axis(counter_2, dimension, coordinate_axis);
			cross_product_vectors.push_back(coordinate_axis);
		}

		// ~~~~~
		MathVector old_normal_of_hyperface = normal_of_hyperface;
		// ~~~~~

		normal_of_hyperface = MathVector::crossProduct(cross_product_vectors);

		// ~~~~~
		double eps = 1e-5;
		if (normal_of_hyperface * normal_of_hyperface < eps)
			normal_of_hyperface = old_normal_of_hyperface;
		else
			normal_of_hyperface.normalize();
		// ~~~~~
	}
}

// all_points - копирование, поскольку они будут изменяться в несимплициальном случае (было const ...&)
void wrapping_algorithm(std::vector<MathVector> all_points, std::unordered_set<std::unordered_set<size_t>, decltype(&us_hash)>& convex_hull)
{
	// -----
	convex_hull.clear();
	// -----

	// Initialize remaining indexes of points, which can be considered.
	std::unordered_set<size_t> interest_indexes_of_points;
	for (size_t counter = 0; counter < all_points.size(); ++counter)
		interest_indexes_of_points.insert(counter);

	// Declare queue of tagged hyperfaces and set of bounding subfaces.
	std::queue<std::unordered_set<size_t>> queue_of_hyperfaces;
	std::unordered_set<std::unordered_set<size_t>, decltype(&us_hash)> bounding_subfaces(0, us_hash);

	// Find first hyperface.
	std::unordered_set<size_t> current_hyperface;
	
	// !!!!!
	find_first_hyperface(all_points, current_hyperface);
	//current_hyperface = std::unordered_set<size_t>({ 0, 2, 6 });
	// !!!!!

	// Push first hyperface in queue and find its subfaces.
	queue_of_hyperfaces.push(current_hyperface);
	find_subfaces(current_hyperface, bounding_subfaces);

	std::unordered_set<std::unordered_set<size_t>, decltype(&us_hash)> subfaces_of_current_hyperface(0, us_hash);

	// -----
	std::unordered_set<std::unordered_set<size_t>, decltype(&us_hash)> subfaces(0, us_hash);
	std::unordered_set<size_t> new_hyperface;
	size_t new_vertex_index;
	// -----	

	std::unordered_set<std::unordered_set<size_t>, decltype(&us_hash)> intersection_of_subfaces(0, us_hash);
	std::unordered_set<std::unordered_set<size_t>, decltype(&us_hash)> subfaces_of_new_hyperface(0, us_hash);

	std::vector<MathVector> cross_product_vectors;

	MathVector normal_of_hyperface;
	MathVector normal_of_subface;

	while (!queue_of_hyperfaces.empty())
	{
		// Pop hyperface from queue and find its subfaces.
		current_hyperface = queue_of_hyperfaces.front();
		queue_of_hyperfaces.pop();
		find_subfaces(current_hyperface, subfaces_of_current_hyperface);

		// -----
		// Find intersection of subfaces_of_current_hyperface and bounding_subfaces.
		// ?????
		subfaces.clear();
		// ?????
		for (auto subface : subfaces_of_current_hyperface)
			if (bounding_subfaces.count(subface))
				subfaces.insert(subface);
		// -----
		
		// Find n.
		cross_product_vectors.clear();

		for (auto it = std::next(current_hyperface.begin(), 1); it != current_hyperface.end(); ++it)
			cross_product_vectors.push_back(all_points[*it] - all_points[*current_hyperface.begin()]);

		normal_of_hyperface = MathVector::crossProduct(cross_product_vectors);
		normal_of_hyperface.normalize();

		// @@@@@
		// Ориентируем нормаль гиперграни в полупространство под aff(F)
		// Это не подходит для вырожденного случая (несимплициального), т.к. там ск. пр-ие может быть == 0
		for (size_t counter = 0; counter < all_points.size(); ++counter)
			if (!current_hyperface.count(counter))
			{
				if ((all_points[counter] - all_points[*current_hyperface.begin()]) * normal_of_hyperface < 0)
					normal_of_hyperface = -1.0 * normal_of_hyperface;
				break;
			}

		// @@@@@


		for (auto subface : subfaces)
		{
			// Find a.
			cross_product_vectors.clear();
			// -----
			cross_product_vectors.push_back(normal_of_hyperface);
			// -----
			for (auto it = std::next(subface.begin(), 1); it != subface.end(); ++it)
				cross_product_vectors.push_back(all_points[*it] - all_points[*subface.begin()]);

			normal_of_subface = MathVector::crossProduct(cross_product_vectors);
			normal_of_subface.normalize();

			// @@@@@
			// (Возможно можно как-то неявно задать ориентацию вектора a)
			// Ориентируем нормаль подграни в противоположном направлении от вершины данной гиперграни вне текущей подграни
			for (auto index : current_hyperface)
				if (!subface.count(index))
				{
					if ((all_points[index] - all_points[*subface.begin()]) * normal_of_subface > 0)
						normal_of_subface = -1.0 * normal_of_subface;
					break;
				}

			
			// @@@@@


			// !!!!!



			// Find vertex of a new hyperface.
			// ?????
			interest_indexes_of_points.clear();
			for (size_t counter = 0; counter < all_points.size(); ++counter)
				interest_indexes_of_points.insert(counter);

			for (auto vertex : current_hyperface)
				interest_indexes_of_points.erase(vertex);
			// ?????
			new_vertex_index = wrapping(false, normal_of_hyperface, all_points[*subface.begin()], normal_of_subface, interest_indexes_of_points, all_points);
			
			// Add new hyperface in queue.
			new_hyperface = subface;
			new_hyperface.insert(new_vertex_index);
			queue_of_hyperfaces.push(new_hyperface);
			// -----

			// Recalculate set of bounding subfaces
			// ~~~~~
			// 1) intersection_of_subfaces = Г П subfaces(F')
			// Г = bounding_subfaces; subfaces(F') = subfaces_of_new_hyperface

			find_subfaces(new_hyperface, subfaces_of_new_hyperface);
			
			// !!!!!
			intersection_of_subfaces.clear();
			for (auto subface_of_new_hyperface : subfaces_of_new_hyperface)
				if (bounding_subfaces.count(subface_of_new_hyperface))
					intersection_of_subfaces.insert(subface_of_new_hyperface);
			// !!!!!

			// 2) Г = Г U F'
			for (auto subface_of_new_hyperface : subfaces_of_new_hyperface)
				bounding_subfaces.insert(subface_of_new_hyperface);

			// 3) Г = Г \ intersection_of_subfaces
			for (auto subface_from_intersection : intersection_of_subfaces)
				bounding_subfaces.erase(subface_from_intersection);
			// ~~~~~

			
		}

		//convex_hull.insert(new_hyperface);
		convex_hull.insert(current_hyperface);
	}
}

// -----

int signOfSemiSpace(const std::vector<MathVector>& points_of_hyperface, const MathVector& test_point)
{
	std::vector<MathVector> cross_product_vectors;
	for (auto it = std::next(points_of_hyperface.begin(), 1); it != points_of_hyperface.end(); ++it)
		cross_product_vectors.push_back(*it - *points_of_hyperface.begin());

	MathVector coefficients = MathVector::crossProduct(cross_product_vectors);

	double res = 0.0;
	for (size_t counter = 0; counter < test_point.getDimension(); ++counter)
		res += coefficients[counter] * (test_point[counter] - (*points_of_hyperface.begin())[counter]);

	return res > 0 ? 1 : res < 0 ? -1 : 0;
}

// ~~~~~
void testPolyhedron(size_t num_of_interior_points, const std::vector<MathVector>& vertices, const std::vector<std::vector<MathVector>>& hyperfaces)
{
	std::vector<MathVector> test_points = vertices;
	
	// !!!!!
	size_t dimension = test_points.begin()->getDimension();
	// !!!!!

	// -----
	// Find MIN and MAX of each coordinate.

	std::vector<double> min_coordinate(dimension);
	std::vector<double> max_coordinate(dimension);

	for (size_t coordinate = 0; coordinate < dimension; ++coordinate)
	{
		min_coordinate[coordinate] = max_coordinate[coordinate] = (*test_points.begin())[coordinate];

		for (auto point : test_points)
		{
			if (point[coordinate] < min_coordinate[coordinate])
				min_coordinate[coordinate] = point[coordinate];
			else
				if (point[coordinate] > max_coordinate[coordinate])
					max_coordinate[coordinate] = point[coordinate];
		}
	}
			
	// -----

	boost::mt19937 generator;
	std::vector<boost::uniform_real<>> distributions(dimension);
	for (size_t coordinate = 0; coordinate < dimension; ++coordinate)
		distributions[coordinate] = boost::uniform_real<>(min_coordinate[coordinate], max_coordinate[coordinate]);

	MathVector interior_point(dimension);
	std::vector<int> signs;

	for (size_t counter_1 = 0; counter_1 < num_of_interior_points; ++counter_1)
	{
		while (true)
		{
			for (size_t coordinate = 0; coordinate < dimension; ++coordinate)
				interior_point[coordinate] = distributions[coordinate](generator);


			for (auto hyperface : hyperfaces)
				signs.push_back(signOfSemiSpace(hyperface, interior_point));

			// -----
			// Мне не нравится такой подход (поискать в будущем что-то более высокоуровневое из std::algorithm (а-ля алгоритмы stl))
			bool isEqualSigns = true;

			for (auto sign : signs)
				if (sign * signs[0] < 0)
				{
					isEqualSigns = false;
					break;
				}

			signs.clear();
			// -----

			if (isEqualSigns)
			{
				// =====
				std::cout << interior_point << std::endl;
				// =====

				test_points.push_back(interior_point);
				break;
			}

		}
	}

	std::unordered_set<std::unordered_set<size_t>, decltype(&us_hash)> convex_hull(0, us_hash);
	wrapping_algorithm(test_points, convex_hull);

	// -----
	//std::ofstream out_points("C:\\Users\\Алексей\\PycharmProjects\\test\\data\\points.txt");
	std::ofstream out_points("data\\points.txt");

	for (auto point : test_points)
		out_points << point << std::endl;

	out_points.close();


	//std::ofstream out_hyperfaces("C:\\Users\\Алексей\\PycharmProjects\\test\\data\\faces.txt");
	std::ofstream out_hyperfaces("data\\faces.txt");

	for (auto hyperface : convex_hull)
	{
		for (auto index : hyperface)
			out_hyperfaces << index << " ";

		out_hyperfaces << std::endl;
	}

	out_hyperfaces.close();

	// -----

	// $$$$$
	for (auto hyperface : convex_hull)
	{
		std::cout << "{ ";
		for (auto index : hyperface)
			std::cout << index << " ";
		std::cout << "}" << std::endl;
	}
	// $$$$$

	system("python plot_3D4D.py");
}

void testRandomPointCloud(size_t num_of_points)
{

}

int main()
{
	// -----
	size_t num_of_interior_points = 0;

	// **** // 
	//  3D  //
	// **** //
	
	// ~~~~~~~
	// Pyramid
	// ~~~~~~~
	/*
	std::vector<MathVector> vertices;
	vertices.push_back(MathVector({ 0.0, 0.0, 1.0 }));
	vertices.push_back(MathVector({ -0.5, pow(3, 0.5) / 2.0, 0.0 }));
	vertices.push_back(MathVector({ -0.5, -pow(3, 0.5) / 2.0, 0.0 }));
	vertices.push_back(MathVector({ 1.0, 0.0, 0.0 }));

	std::vector<std::vector<MathVector>> hyperfaces;
	hyperfaces.push_back(std::vector<MathVector>({ vertices[0], vertices[1], vertices[2] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[0], vertices[3], vertices[1] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[0], vertices[2], vertices[3] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[2], vertices[1], vertices[3] }));
	*/
	
	// +++++

	// ~~~~~~~~~~
	// Octahedron
	// ~~~~~~~~~~
	/*
	std::vector<MathVector> vertices;
	vertices.push_back(MathVector({ 1.0, 0.0, 0.0 }));
	vertices.push_back(MathVector({ -1.0, 0.0, 0.0 }));
	vertices.push_back(MathVector({ 0.0, 1.0, 0.0 }));
	vertices.push_back(MathVector({ 0.0, -1.0, 0.0 }));
	vertices.push_back(MathVector({ 0.0, 0.0, 1.0 }));
	vertices.push_back(MathVector({ 0.0, 0.0, -1.0 }));

	std::vector<std::vector<MathVector>> hyperfaces;
	hyperfaces.push_back(std::vector<MathVector>({ vertices[1], vertices[4], vertices[2] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[2], vertices[4], vertices[0] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[2], vertices[0], vertices[5] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[1], vertices[2], vertices[5] }));

	hyperfaces.push_back(std::vector<MathVector>({ vertices[4], vertices[1], vertices[3] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[4], vertices[3], vertices[0] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[1], vertices[5], vertices[3] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[3], vertices[5], vertices[0] }));
	*/
	

	// +++++

	// ~~~~~~~~~~~
	// Icosahedron
	// ~~~~~~~~~~~
	/*
	std::vector<MathVector> vertices;

	vertices.push_back(MathVector({ -0.692, 0.000, 0.427 }));
	vertices.push_back(MathVector({ 0.000, 0.427, -0.692 }));
	vertices.push_back(MathVector({ 0.000, 0.427, 0.692 }));
	vertices.push_back(MathVector({ 0.692, 0.000, -0.427 }));
	vertices.push_back(MathVector({ -0.427, -0.692, 0.000 }));
	vertices.push_back(MathVector({ -0.427, 0.692, 0.000 }));
	vertices.push_back(MathVector({ 0.000, -0.427, 0.692 }));
	vertices.push_back(MathVector({ 0.427, 0.692, 0.000 }));
	vertices.push_back(MathVector({ 0.000, -0.427, -0.692 }));
	vertices.push_back(MathVector({ 0.692, 0.000, 0.427 }));
	vertices.push_back(MathVector({ 0.427, -0.692, 0.000 }));
	vertices.push_back(MathVector({ -0.692, 0.000, -0.427 }));

	std::vector<std::vector<MathVector>> hyperfaces;
	
	hyperfaces.push_back(std::vector<MathVector>({ vertices[9], vertices[2], vertices[6] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[1], vertices[11], vertices[5] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[11], vertices[1], vertices[8] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[0], vertices[11], vertices[4] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[3], vertices[1], vertices[7] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[3], vertices[8], vertices[1] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[9], vertices[3], vertices[7] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[0], vertices[6], vertices[2] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[4], vertices[10], vertices[6] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[1], vertices[5], vertices[7] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[7], vertices[5], vertices[2] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[8], vertices[3], vertices[10] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[4], vertices[11], vertices[8] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[9], vertices[7], vertices[2] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[10], vertices[9], vertices[6] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[0], vertices[5], vertices[11] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[0], vertices[2], vertices[5] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[8], vertices[10], vertices[4] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[3], vertices[9], vertices[10] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[6], vertices[0], vertices[4] }));
	*/
	// -----

	// ~~~~
	// Cube
	// ~~~~
	/*
	std::vector<MathVector> vertices;

	vertices.push_back(MathVector({ -1.0, -1.0, -1.0 }));
	vertices.push_back(MathVector({ -1.0,  1.0, -1.0 }));
	vertices.push_back(MathVector({  1.0, -1.0, -1.0 }));
	vertices.push_back(MathVector({  1.0,  1.0, -1.0 }));
	vertices.push_back(MathVector({ -1.0, -1.0,  1.0 }));
	vertices.push_back(MathVector({ -1.0,  1.0,  1.0 }));
	vertices.push_back(MathVector({  1.0, -1.0,  1.0 }));
	vertices.push_back(MathVector({  1.0,  1.0,  1.0 }));

	std::vector<std::vector<MathVector>> hyperfaces;
	
	
	hyperfaces.push_back(std::vector<MathVector>({ vertices[6], vertices[2], vertices[7] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[3], vertices[7], vertices[2] }));

	hyperfaces.push_back(std::vector<MathVector>({ vertices[7], vertices[3], vertices[5] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[1], vertices[5], vertices[3] }));
	
	hyperfaces.push_back(std::vector<MathVector>({ vertices[5], vertices[1], vertices[4] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[0], vertices[4], vertices[1] }));

	hyperfaces.push_back(std::vector<MathVector>({ vertices[4], vertices[0], vertices[6] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[2], vertices[6], vertices[0] }));

	hyperfaces.push_back(std::vector<MathVector>({ vertices[4], vertices[6], vertices[5] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[7], vertices[5], vertices[6] }));

	hyperfaces.push_back(std::vector<MathVector>({ vertices[0], vertices[1], vertices[2] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[3], vertices[2], vertices[1] }));
	*/
	
	/*
	hyperfaces.push_back(std::vector<MathVector>({ vertices[0], vertices[1], vertices[3], vertices[2] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[4], vertices[6], vertices[7], vertices[5] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[6], vertices[2], vertices[3], vertices[7] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[7], vertices[3], vertices[1], vertices[5] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[5], vertices[1], vertices[4], vertices[0] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[4], vertices[0], vertices[2], vertices[6] }));
	*/

	// **** // 
	//  4D  //
	// **** //

	// ~~~~~~
	// 5-cell
	// ~~~~~~

/*
	std::vector<MathVector> vertices;

	vertices.push_back(MathVector({  1.0,  1.0,  1.0, 0.0 }));
	vertices.push_back(MathVector({  1.0, -1.0, -1.0, 0.0 }));
	vertices.push_back(MathVector({ -1.0,  1.0, -1.0, 0.0 }));
	vertices.push_back(MathVector({ -1.0, -1.0,  1.0, 0.0 }));
	vertices.push_back(MathVector({  0.0,  0.0,  0.0, sqrt(5.0) }));

	
	std::vector<std::vector<MathVector>> hyperfaces;

	hyperfaces.push_back(std::vector<MathVector>({ vertices[1], vertices[0], vertices[2], vertices[3] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[4], vertices[0], vertices[2], vertices[1] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[4], vertices[0], vertices[1], vertices[3] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[4], vertices[0], vertices[3], vertices[2] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[4], vertices[1], vertices[2], vertices[3] }));
*/
	// -----
	
	// ~~~~~~~
	// 16-cell
	// ~~~~~~~
/*
	std::vector<MathVector> vertices;

	vertices.push_back(MathVector({  1.0,  0.0,  0.0, 0.0 }));
	vertices.push_back(MathVector({ -1.0,  0.0,  0.0, 0.0 }));

	vertices.push_back(MathVector({  0.0,  1.0,  0.0, 0.0 }));
	vertices.push_back(MathVector({  0.0, -1.0,  0.0, 0.0 }));

	vertices.push_back(MathVector({  0.0,  0.0,  1.0, 0.0 }));
	vertices.push_back(MathVector({  0.0,  0.0, -1.0, 0.0 }));

	vertices.push_back(MathVector({ 0.0,  0.0,  0.0,  1.0 }));
	vertices.push_back(MathVector({ 0.0,  0.0,  0.0, -1.0 }));


	std::vector<std::vector<MathVector>> hyperfaces;

	
	hyperfaces.push_back(std::vector<MathVector>({ vertices[6], vertices[5], vertices[1], vertices[2] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[6], vertices[5], vertices[2], vertices[0] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[6], vertices[5], vertices[0], vertices[3] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[6], vertices[5], vertices[3], vertices[1] }));
	
	hyperfaces.push_back(std::vector<MathVector>({ vertices[6], vertices[4], vertices[2], vertices[1] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[6], vertices[4], vertices[0], vertices[2] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[6], vertices[4], vertices[3], vertices[0] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[6], vertices[4], vertices[1], vertices[3] }));


	hyperfaces.push_back(std::vector<MathVector>({ vertices[5], vertices[7], vertices[1], vertices[2] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[5], vertices[7], vertices[2], vertices[0] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[5], vertices[7], vertices[0], vertices[3] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[5], vertices[7], vertices[3], vertices[1] }));

	hyperfaces.push_back(std::vector<MathVector>({ vertices[4], vertices[7], vertices[2], vertices[1] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[4], vertices[7], vertices[0], vertices[2] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[4], vertices[7], vertices[3], vertices[0] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[4], vertices[7], vertices[1], vertices[3] }));
	// -----
*/
	// ~~~~~~~~~
	// Tesseract
	// ~~~~~~~~~

	std::vector<MathVector> vertices;

	vertices.push_back(MathVector({ -1.0, -1.0, -1.0, -1.0 }));
	vertices.push_back(MathVector({ -1.0, -1.0, -1.0,  1.0 }));
	vertices.push_back(MathVector({ -1.0, -1.0,  1.0, -1.0 }));
	vertices.push_back(MathVector({ -1.0, -1.0,  1.0,  1.0 }));
	vertices.push_back(MathVector({ -1.0,  1.0, -1.0, -1.0 }));
	vertices.push_back(MathVector({ -1.0,  1.0, -1.0,  1.0 }));
	vertices.push_back(MathVector({ -1.0,  1.0,  1.0, -1.0 }));
	vertices.push_back(MathVector({ -1.0,  1.0,  1.0,  1.0 }));
	vertices.push_back(MathVector({  1.0, -1.0, -1.0, -1.0 }));
	vertices.push_back(MathVector({  1.0, -1.0, -1.0,  1.0 }));
	vertices.push_back(MathVector({  1.0, -1.0,  1.0, -1.0 }));
	vertices.push_back(MathVector({  1.0, -1.0,  1.0,  1.0 }));
	vertices.push_back(MathVector({  1.0,  1.0, -1.0, -1.0 }));
	vertices.push_back(MathVector({  1.0,  1.0, -1.0,  1.0 }));
	vertices.push_back(MathVector({  1.0,  1.0,  1.0, -1.0 }));
	vertices.push_back(MathVector({  1.0,  1.0,  1.0,  1.0 }));

	std::vector<std::vector<MathVector>> hyperfaces;

	// Inner cube
	hyperfaces.push_back(std::vector<MathVector>({ vertices[10], vertices[8], vertices[12], vertices[14] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[14], vertices[12], vertices[4], vertices[6] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[6], vertices[4], vertices[0], vertices[2] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[2], vertices[0], vertices[8], vertices[10] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[2], vertices[10], vertices[14], vertices[6] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[8], vertices[0], vertices[4], vertices[12] }));

	// Outer cube
	hyperfaces.push_back(std::vector<MathVector>({ vertices[11], vertices[9], vertices[13], vertices[15] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[15], vertices[13], vertices[5], vertices[7] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[7], vertices[5], vertices[1], vertices[3] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[3], vertices[1], vertices[9], vertices[11] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[3], vertices[11], vertices[15], vertices[7] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[9], vertices[1], vertices[5], vertices[13] }));

	// Bottom partitions
	hyperfaces.push_back(std::vector<MathVector>({ vertices[8], vertices[9], vertices[13], vertices[12] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[12], vertices[13], vertices[5], vertices[4] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[4], vertices[5], vertices[1], vertices[0] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[0], vertices[1], vertices[9], vertices[8] }));

	// Upper partitions
	hyperfaces.push_back(std::vector<MathVector>({ vertices[10], vertices[11], vertices[15], vertices[14] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[14], vertices[15], vertices[7], vertices[6] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[6], vertices[7], vertices[3], vertices[2] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[2], vertices[3], vertices[11], vertices[10] }));

	// Side partitions
	hyperfaces.push_back(std::vector<MathVector>({ vertices[11], vertices[10], vertices[8], vertices[9] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[15], vertices[14], vertices[12], vertices[13] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[7], vertices[6], vertices[4], vertices[5] }));
	hyperfaces.push_back(std::vector<MathVector>({ vertices[3], vertices[2], vertices[0], vertices[1] }));


	testPolyhedron(num_of_interior_points, vertices, hyperfaces);

	system("pause");
	return 0;
}