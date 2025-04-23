
// MIT License
//
// Copyright (c) 2019 SuYuxi
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// Comment by Gregor Skok:
// The original code by SuYuxi has beed modified a lot to make it faster and to include new features

#include <vector>
#include <queue>
#include <memory>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace kdtree
	{
	using namespace std;

	typedef float PointType;
	//typedef vector<PointType> Point; //Presents Point type by vector<int> like (x, y, z)

	struct Point_str  // Point_str Structure declaration
		{
		vector<PointType> coords;
		size_t index;

		void set(const vector<PointType> &coords_, const size_t &index_)
			{
			coords = coords_;
			index = index_;
			}

		void free_memory()
			{
			free_vector_memory(coords);
			}

		string output() const
			{
			ostringstream s1;
			s1 << "(";
			for (unsigned long il=0; il < coords.size(); il++)
				{
				s1 << coords[il];
				if (il< coords.size() - 1)
					s1 << ",";
				}
			s1 << "):" << index;
			return(s1.str());
			}
 		};


	bool operator==(const Point_str& lhs, const Point_str& rhs)
		{
		if (lhs.coords == rhs.coords && lhs.index == rhs.index)
			return(true);
		return(false);
		}

	bool operator!=(const Point_str& lhs, const Point_str& rhs)
		{
		if (lhs.coords != rhs.coords || lhs.index != rhs.index)
			return(true);
		return(false);
		}

	struct KdTreeNode
		{
		KdTreeNode(){};

		KdTreeNode(const Point_str &_val)
			{
			val = _val;
			leftNode = nullptr;
			rightNode = nullptr;
			}

		Point_str val;
		vector<PointType> coords_max;
		vector<PointType> coords_min;
		KdTreeNode * leftNode;
		KdTreeNode * rightNode;

		void free_memory()
			{
			val.free_memory();
			free_vector_memory(coords_max);
			free_vector_memory(coords_min);
			delete leftNode;
			leftNode = nullptr;
			delete rightNode;
			rightNode = nullptr;
			}

		string output() const
			{
			ostringstream s1;
			s1 << "(";
			for (unsigned long il=0; il < val.coords.size(); il++)
				{
				s1 << val.coords[il];
				if (il< val.coords.size() - 1)
					s1 << ",";
				}
			s1 << "):" << val.index ;
			s1 << " min(" << output_vector_as_string(coords_min, ",") << ") max(" <<  output_vector_as_string(coords_max, ",") << ")";
			return(s1.str());
			}

		};

	class KdTree
		{
		public:
		int32_t dimension;
		KdTreeNode * root;

		KdTree()
			{
			dimension = 0;
			root = nullptr;
			}

		KdTree(vector<Point_str> &points)
			{
			KdTree();
			buildKdTree(points);
			}

		~KdTree()
			{
			free_memory();
			}

		bool is_empty() const
			{
			if (root == nullptr)
				return true;
			return false;
			}

		bool buildKdTree(vector<Point_str> &points)
			{ //All the points must be of same dimension like (x, y ,z) or (x, y) or what ever you like
			if(points.empty() | points[0].coords.empty()) return false;
			dimension = points[0].coords.size();
			root = buildHelper(points, 0, points.size() - 1, 0);


			// add Bounding Box data for faster NN search
			add_Bounding_Box_data_to_the_tree();

			return root != nullptr ? true : false; //return false if cannot build kd Tree by input points;
			}


		KdTreeNode * buildHelper(vector<Point_str>& points, const int32_t leftBorder, const int32_t rightBorder, const int32_t depth)
			{ //recursively create kd Tree
			if(leftBorder > rightBorder) return nullptr;
			int32_t curDim = depth % dimension;
			int32_t midInx = leftBorder + (rightBorder - leftBorder) / 2;
			//sort(points.begin() + leftBorder, points.begin() + rightBorder + 1, [curDim](const Point_str& a, const Point_str& b) { return a[curDim] < b[curDim]; });

			// ----------------------
			// ni nujno uporabiti celotnega sorta - namesto tega se lahko uporabi nth_element z še nekaj dodatnega dela - to je potem kar hitrejse - upam da zadeva dela prav
			nth_element(points.begin() + leftBorder, points.begin() + midInx,  points.begin() + rightBorder + 1, [curDim](const Point_str& a, const Point_str& b) { return a.coords[curDim] < b.coords[curDim]; });
			// sedaj je treba še na levi strani skupaj prestaviti vrednosti, ki so enake kot midInx (nth_element funkcija tega ne garantira)
			if (midInx - leftBorder > 1)
				{
				int32_t last_taken_indx=midInx;
				if (points[midInx-1].coords[curDim] == points[midInx].coords[curDim])
					last_taken_indx--;

				for (long il=midInx-2; il >= leftBorder; il--)
					if (points[il].coords[curDim] == points[midInx].coords[curDim])
						{
						swap(points[il],points[last_taken_indx-1]);
						last_taken_indx--;
						}
				}
			// ----------------------



			while(midInx > leftBorder && points[midInx - 1].coords[curDim] == points[midInx].coords[curDim]) // keep all the points with point[splitDim] >= midInx[splitDim] on the right of midInx
			{
				midInx -= 1;
			}

			/*cout << "depth: " << depth << endl;
			cout << "curDim: " << curDim << endl;
			for (int32_t il=leftBorder; il <= rightBorder; il++)
				{
				cout << points[il].output();
				if (il == midInx) cout << "*";
				cout << endl;
				}
			cout << "------" << endl;
			*/
			KdTreeNode * node = new KdTreeNode(points[midInx]);

			/*
			// calculate max and min coordinate values in the leaf for the Bounding Box data
			vector <PointType> min_coords = points[leftBorder].coords;
			vector <PointType> max_coords = points[leftBorder].coords;
			for (int32_t il=leftBorder + 1; il <= rightBorder; il++)
			for (int32_t ic=0; ic < dimension; ic ++)
				{
				if (points[il].coords[ic] < min_coords[ic]) min_coords[ic] = points[il].coords[ic];
				if (points[il].coords[ic] > max_coords[ic]) max_coords[ic] = points[il].coords[ic];
				}
			node->coords_max=max_coords;
			node->coords_min=min_coords;
			*/
			node->leftNode = buildHelper(points, leftBorder, midInx - 1, depth + 1);
			node->rightNode = buildHelper(points, midInx + 1, rightBorder, depth + 1);

			return node;
			}

		void add_Bounding_Box_data_to_the_tree()
			{
			vector <PointType> coords_min (dimension, 0);
			vector <PointType> coords_max (dimension, 0);
			add_Bounding_Box_data_to_the_tree_helper(root, coords_min, coords_max);
			}

		void add_Bounding_Box_data_to_the_tree_helper(KdTreeNode * const node, vector <PointType> &coords_min_out, vector <PointType> &coords_max_out)
			{
			node->coords_min=node->val.coords;
			node->coords_max=node->val.coords;

			//cout << node->output() << endl;

			if (node->leftNode != nullptr)
				add_Bounding_Box_data_to_the_tree_helper(node->leftNode, node->coords_min, node->coords_max);
			if (node->rightNode != nullptr)
				add_Bounding_Box_data_to_the_tree_helper(node->rightNode, node->coords_min, node->coords_max);

			for (int32_t ic=0; ic < dimension; ic ++)
				{
				coords_min_out[ic]=min(coords_min_out[ic],node->coords_min[ic]);
				coords_max_out[ic]=max(coords_max_out[ic],node->coords_max[ic]);
				}

			//cout << "-" << node->output() << endl;
			}

		bool buildKdTree_and_do_not_change_the_kdtree_points_vector(const vector<Point_str> &kdtree_points)
			{
			vector <kdtree::Point_str> kdtree_points_temp = kdtree_points;
			bool result = buildKdTree(kdtree_points_temp);
			free_vector_memory(kdtree_points_temp); // free memory
			return(result);
			}

		void free_memory()
			{
			if (root != nullptr)
				free_memory_helper(root);
			dimension = 0;
			}

		void free_memory_helper(KdTreeNode * &node)
			{
			if (node->leftNode != nullptr)
				free_memory_helper(node->leftNode);
			if (node->rightNode != nullptr)
				free_memory_helper(node->rightNode);
			delete node;
			node = nullptr;
			}

		void printKdTree2() const
			{
			if (root != nullptr)
				printKdTree2_helper(root, 0, 0, true);
			else
				cout << "Kdtree does not contain any nodes ! " << endl;
			}

		void printKdTree2_helper(const KdTreeNode * const node, const int32_t depth, const int indent, const bool is_left) const
			{
			int delta_indent=12;
			for (int32_t il=0; il < (depth-1)*delta_indent; il++)
				cout << " ";
			if (depth > 0)
				for (int32_t il=0; il < delta_indent; il++)
					cout << "-";

			int32_t curDim = depth % dimension;

			cout << "L" << depth;
			if (depth > 0)
				{
				if (is_left) cout << "l";
				else cout << "r";
				}
			cout << ":";

			cout << "(";
			for (unsigned long il=0; il < node->val.coords.size(); il++)
				{
				if ((int32_t)il == curDim)
					cout << "*";
				cout << node->val.coords[il];
				if ( il < node->val.coords.size() - 1)
					cout << ",";
				}
			cout << "):" << node->val.index << " min(" << output_vector_as_string(node->coords_min, ",") << ") max(" <<  output_vector_as_string(node->coords_max, ",") << ")"<< endl;


			if(node->leftNode != nullptr)
				printKdTree2_helper(node->leftNode, depth +1, indent + delta_indent, true);
			else
				{
				for (int32_t il=0; il < (depth)*delta_indent; il++)
					cout << " ";
				for (int32_t il=0; il < delta_indent; il++)
					cout << "-";
				cout << "L" << depth+1 << ":""nullptr" << endl;
				}

			if(node->rightNode != nullptr)
				printKdTree2_helper(node->rightNode, depth +1, indent + delta_indent, false);
			else
				{
				for (int32_t il=0; il < (depth)*delta_indent; il++)
					cout << " ";
				for (int32_t il=0; il < delta_indent; il++)
					cout << "-";
				cout << "L" << depth+1 << ":""nullptr" << endl;
				}
			//indent = indent_old;
			}


		//print all Kd Tree's nodes layer by layer using breadth first search
		size_t count_number_of_nonnull_nodes() const
			{
			size_t counter=0;
			queue<KdTreeNode *> q;
			q.emplace(root);
			KdTreeNode * node;
			while(!q.empty())
				{
					node = q.front();
					q.pop();
					if(node != nullptr)
					{
						counter++;
						//cout << node->val.output() << endl;
						//cout << "(";
						//for_each(node->val.coords.begin(), node->val.coords.end() - 1, [](PointType& num) { cout << num << ", ";});
						//cout << *(node->val.coords.end() - 1) << ")" << endl;
						q.emplace(node->leftNode);
						q.emplace(node->rightNode);
					}
				}
			return(counter);
			}


		// rebuild a balanced kd-tree from scratch using all the points in the current tree - delanje tega med racunanjem PADa ni pohitrilo zadeve (2024-04)
		bool rebuild_a_balanced_KdTree()
			{
			vector<Point_str> points;
			// get all points
			rebuild_a_balanced_KdTree_helper(root, points);
			// rebuild the tree
			return(buildKdTree(points));
			}

		void rebuild_a_balanced_KdTree_helper( KdTreeNode * &node, vector<Point_str> &points)
			{
			ERRORIF(node == nullptr);

			if(node->leftNode != nullptr)
				rebuild_a_balanced_KdTree_helper(node->leftNode, points);
			if(node->rightNode != nullptr)
				rebuild_a_balanced_KdTree_helper(node->rightNode, points);
			// save this point
			points.push_back(node->val);
			// free memory
			delete node;
			node = nullptr;
			}


		const KdTreeNode * findMin(const KdTreeNode * const node, const int32_t dim, const int32_t depth) const //find the node with the minimum value on dim dimension from depth
		{
			if(root == nullptr) return nullptr;
			const KdTreeNode * minimum = node;
			findMinHelper(node, minimum, dim, depth);
			return minimum;
		}

		void findMinHelper(const KdTreeNode * const node, const KdTreeNode * &minimum, const int32_t dim, const int32_t depth) const
			{
			if(node == nullptr) return;
			if(node->val.coords[dim] < minimum->val.coords[dim]) { minimum = node; }
			findMinHelper(node->leftNode, minimum, dim, depth + 1);
			if(depth % dimension != dim) { findMinHelper(node->rightNode, minimum, dim, depth + 1); }
			}

		// Tree will not be balanced
		bool addNode(const Point_str & point)
			{

			// special case if tree is empty
			if (root == nullptr)
				{
				vector <Point_str> points = {point};
				if (buildKdTree(points))
					return true;
				return false;
				}

			if(point.coords.size() != (unsigned long)dimension) return false;

			KdTreeNode * node = root;

			int32_t depth = 0;
			int32_t curDim;
			while(true)
				{
				curDim = depth % dimension;
				// update BB data
				for (int32_t ic=0; ic < dimension; ic ++)
					{
					node->coords_min[ic]=min(node->coords_min[ic],point.coords[ic]);
					node->coords_max[ic]=max(node->coords_max[ic],point.coords[ic]);
					}

				if(point.coords[curDim] < node->val.coords[curDim])
					{
					if(node->leftNode == nullptr)
						{
						node->leftNode = new KdTreeNode(point);
						// add BB data
						node->leftNode->coords_min = point.coords;
						node->leftNode->coords_max = point.coords;
						return true;
						}
					node = node->leftNode;
					}
				else
					{
					if(node->rightNode == nullptr)
						{
						node->rightNode = new KdTreeNode(point);
						// add BB data
						node->rightNode->coords_min = point.coords;
						node->rightNode->coords_max = point.coords;
						return true;
						}
					node = node->rightNode;
					}
				depth += 1;
				}

			return false;
			}

		void deleteNode(const Point_str &point)
			{
			if(root == nullptr || point.coords.size() != (unsigned long)dimension) return;
			if(deleteNodeHelper(root, point, 0))
				{
				delete root;
				root = nullptr;
				}
			}

		bool deleteNodeHelper(KdTreeNode * const node, const Point_str& point, const int32_t depth)
			{ //return true means the child node should be deleted
			if(node == nullptr) return false;
			int32_t curDim = depth % dimension;
			if(point == node->val)
				{
				if(node->rightNode != nullptr)
					{
					const KdTreeNode * minimumNode = findMin(node->rightNode, curDim, depth + 1);
					node->val = minimumNode->val; // do not swap(node->val, minimumNode) which would break the structure of kd-tree and result in not finding the node to delete
					if(deleteNodeHelper(node->rightNode, minimumNode->val, depth + 1))
						{
						delete node->rightNode;
						node->rightNode = nullptr;
						}

					// update Bounding Box data
					node->coords_min = node->val.coords;
					node->coords_max = node->val.coords;
					if (node->rightNode != nullptr)
						for (int32_t ic=0; ic < dimension; ic ++)
							{
							node->coords_min[ic]=min(node->rightNode->coords_min[ic],node->coords_min[ic]);
							node->coords_max[ic]=max(node->rightNode->coords_max[ic],node->coords_max[ic]);
							}

					if (node->leftNode != nullptr)
						for (int32_t ic=0; ic < dimension; ic ++)
							{
							node->coords_min[ic]=min(node->leftNode->coords_min[ic],node->coords_min[ic]);
							node->coords_max[ic]=max(node->leftNode->coords_max[ic],node->coords_max[ic]);
							}
					}
				else if(node->leftNode != nullptr)
					{
					const KdTreeNode *  minimumNode = findMin(node->leftNode, curDim, depth + 1);
					node->val = minimumNode->val;
					if(deleteNodeHelper(node->leftNode, minimumNode->val, depth + 1))
						{
						delete node->leftNode;
						node->leftNode = nullptr;
						}

					if (node->rightNode != nullptr)
						{
						delete node->rightNode;
						node->rightNode = nullptr;
						}

					node->rightNode = node->leftNode;
					node->leftNode = nullptr;

					// update Bounding Box data
					node->coords_min = node->val.coords;
					node->coords_max = node->val.coords;
					if (node->rightNode != nullptr)
						for (int32_t ic=0; ic < dimension; ic ++)
							{
							node->coords_min[ic]=min(node->rightNode->coords_min[ic],node->coords_min[ic]);
							node->coords_max[ic]=max(node->rightNode->coords_max[ic],node->coords_max[ic]);
							}
					//cout << "bbb" << endl;
					}
				else
					{
					//cout << "ccc" << endl;
					return true; //return true to inform outter deleteNodeHelper to release this pointer;
					}
				}
			else
				{
				if(point.coords[curDim] < node->val.coords[curDim])
					{
					if(deleteNodeHelper(node->leftNode, point, depth + 1))
						{
						delete node->leftNode;
						node->leftNode = nullptr;
						}
					}
				else
					{
					if(deleteNodeHelper(node->rightNode, point, depth + 1))
						{
						delete node->rightNode;
						node->rightNode = nullptr;
						}
					}
				}
			return false;
			}


		const KdTreeNode * getNode(const Point_str &point) const
			{
			if(root == nullptr || point.coords.size() != (unsigned long)dimension) return nullptr;
			KdTreeNode * node = root;
			int32_t depth = 0;
			int32_t curDim;
			while(node != nullptr)
			{
				if(node->val == point) return node;
				curDim = depth % dimension;
				if(node->val.coords[curDim] > point.coords[curDim])
				{
					node = node->leftNode;
				}
				else
				{
					node = node->rightNode;
				}
				depth += 1;
			}

			return nullptr;
			}


		vector <const KdTreeNode *> get_all_SubNodes(const KdTreeNode * const node) const
			{
			vector <const KdTreeNode *> nodes;
			get_all_SubNodes_Helper(nodes, node);
			return(nodes);
			}

		void get_all_SubNodes_Helper(vector <const KdTreeNode *>  &nodes, const KdTreeNode * const node) const
			{
			if (node != nullptr)
				{
				nodes.push_back(node);
				get_all_SubNodes_Helper(nodes,  node->leftNode);
				get_all_SubNodes_Helper(nodes,  node->rightNode);
				}
			}

		const KdTreeNode * findNearestNode(const Point_str &point) const
			{
			if(root == nullptr || point.coords.size() != (unsigned long)dimension) return nullptr;
			const KdTreeNode * nearestNode = root;
			float minDist_sqr = calDist_sqr(point, nearestNode->val);
			float minDist = sqrt(minDist_sqr);
			findNearestNodeHelper(root, point, minDist_sqr, minDist, nearestNode, 0);
			return nearestNode;
			}

		void findNearestNodeHelper(const KdTreeNode * const node, const Point_str& point, float& minDist_sqr, float& minDist, const KdTreeNode * &nearestNode, const int32_t depth) const
			{
			if(node == nullptr) return;
			int32_t curDim = depth % dimension;
			//float dist = calDist(point, node->val);
			float dist_sqr = calDist_sqr(point, node->val);
			if(dist_sqr < minDist_sqr)
				{
				nearestNode = node;
				minDist_sqr = dist_sqr;
				minDist = sqrt(minDist_sqr);
				}
			//counter_temp++;

			//cout << counter_temp << " " << node->val.output() << " L" << depth << ":" << curDim << " " << minDist << " min(" << output_vector_as_string(node->coords_min, ",") << ") max(" <<  output_vector_as_string(node->coords_max, ",") << ")"<< endl;

			float delta = point.coords[curDim] - node->val.coords[curDim];
			if( delta < 0)
				{
				findNearestNodeHelper(node->leftNode, point, minDist_sqr, minDist, nearestNode, depth + 1);
				if(-delta <= minDist)
					if (test_if_the_hypersphere_and_hyperrectangle_intersect( point.coords, minDist, node->coords_min, node->coords_max ))
						{
						findNearestNodeHelper(node->rightNode, point, minDist_sqr, minDist, nearestNode, depth + 1);
						}
 				}
			else
 				{
				findNearestNodeHelper(node->rightNode, point, minDist_sqr, minDist, nearestNode, depth + 1);
				if(delta < minDist)
					if (test_if_the_hypersphere_and_hyperrectangle_intersect( point.coords, minDist, node->coords_min, node->coords_max ))
						{
						findNearestNodeHelper(node->leftNode, point, minDist_sqr, minDist, nearestNode, depth + 1);
						}
				}

			}

		float calDist_sqr(const Point_str& a, const Point_str& b) const
			{
			float sum = 0;
			for(int32_t inx = 0; inx < dimension; inx++)
				{
				float temp = (float)(b.coords[inx] - a.coords[inx]);
				sum += temp*temp;
				}
			return(sum);
			}

		bool test_if_the_hypersphere_and_hyperrectangle_intersect( const vector <PointType> &circle_center, const float radius, const vector <PointType> &rect_coords_min, const vector <PointType> &rect_coords_max) const
			{
			/*
			// get the point of rectangle that is nearest the circle center
			vector <PointType> nearest_point;
			for (int32_t ic=0; ic < dimension; ic++)
				nearest_point.push_back(clip_value_to_min_max(circle_center[ic],rect_coords_min[ic], rect_coords_max[ic]));

			// the distance between the nearest point and circle center
			float sqr_distance = squared_euclidian_distanceX(circle_center, nearest_point);
			*/

			float sqr_distance = 0;
			for (int32_t ic=0; ic < dimension; ic++)
				{
				float temp = circle_center[ic] - clip_value_to_min_max(circle_center[ic],rect_coords_min[ic], rect_coords_max[ic]);
				sqr_distance+= temp*temp;
				}

			//cout << "- (" <<  output_vector_as_string(nearest_point, ",") << ") " << sqrt(sqr_distance) << endl;

			if (sqr_distance < radius * radius)
				return(true);
			return(false);
			}

		bool test_if_the_hypersphere_and_hyperrectangle_intersect_using_sqr_radius( const vector <PointType> &circle_center, const float sqr_radius, const vector <PointType> &rect_coords_min, const vector <PointType> &rect_coords_max) const
			{
			/*
			// get the point of rectangle that is nearest the circle center
			vector <PointType> nearest_point;
			for (int32_t ic=0; ic < dimension; ic++)
				nearest_point.push_back(clip_value_to_min_max(circle_center[ic],rect_coords_min[ic], rect_coords_max[ic]));

			// the distance between the nearest point and circle center
			float sqr_distance = squared_euclidian_distanceX(circle_center, nearest_point);
			*/

			float sqr_distance = 0;
			for (int32_t ic=0; ic < dimension; ic++)
				{
				float temp = circle_center[ic] - clip_value_to_min_max(circle_center[ic],rect_coords_min[ic], rect_coords_max[ic]);
				/*float temp = 0;
				float cc = circle_center[ic];
				if ( cc < rect_coords_min[ic]) temp= cc - rect_coords_min[ic];
				else if (cc > rect_coords_max[ic]) temp= cc - rect_coords_max[ic];*/
				sqr_distance+= temp*temp;
				}

			//cout << "- (" <<  output_vector_as_string(nearest_point, ",") << ") " << sqrt(sqr_distance) << endl;

			if (sqr_distance < sqr_radius)
				return(true);
			return(false);
			}


		bool test_if_the_hyperrectangle_is_fully_inside_the_sphere( const vector <PointType> &circle_center, const float radius, const vector <PointType> &rect_coords_min, const vector <PointType> &rect_coords_max) const
			{
			float sqr_distance = 0;
			for (int32_t ic=0; ic < dimension; ic++)
				{
				float temp = max(fabs(circle_center[ic] - rect_coords_min[ic]), fabs(circle_center[ic] - rect_coords_max[ic]));
				sqr_distance+= temp*temp;
				}

			if (sqr_distance < radius * radius)
				return(true);
			return(false);
			}

		bool test_if_the_hyperrectangle_is_fully_inside_the_sphere_using_sqr_radius( const vector <PointType> &circle_center, const float sqr_radius, const vector <PointType> &rect_coords_min, const vector <PointType> &rect_coords_max) const
			{
			float sqr_distance = 0;
			for (int32_t ic=0; ic < dimension; ic++)
				{
				float temp = max(fabs(circle_center[ic] - rect_coords_min[ic]), fabs(circle_center[ic] - rect_coords_max[ic]));
				sqr_distance+= temp*temp;
				}

			if (sqr_distance < sqr_radius)
				return(true);
			return(false);
			}


		PointType clip_value_to_min_max (const PointType &val, const PointType &min_limit, const PointType &max_limit) const
			{
			if (val < min_limit) return min_limit;
			if (val > max_limit) return max_limit;
			return(val);
 			}


		vector <const KdTreeNode *> findNearestNodeCluster(const Point_str& point, const float distance) const
			{ //find all nodes with the distance from which to the "point" is less than or equal to the "distance."
			vector <const KdTreeNode *> cluster;
			if(root == nullptr || point.coords.size() != (unsigned long)dimension) return cluster;
			float distance_squared = distance*distance;
			findNearestNodeClusterHelper(root, point, distance, distance_squared, cluster, 0);
			return cluster;
			}


		void findNearestNodeClusterHelper(const KdTreeNode * const node, const Point_str& point, const float distance, const float  distance_sqaured, vector <const KdTreeNode *> &cluster, const int32_t depth) const
			{
			if(node == nullptr) return;

			/*
			if (test_if_the_hyperrectangle_is_fully_inside_the_sphere( point.coords, distance, node->coords_min, node->coords_max))
				{
				vector<std::shared_ptr<KdTreeNode>> subnodes = get_all_SubNodes(node);
				cluster.insert(cluster.end(), subnodes.begin(), subnodes.end());
				}

			else if (test_if_the_hypersphere_and_hyperrectangle_intersect( point.coords, distance, node->coords_min, node->coords_max ))
				{
				float dist_sqr = calDist_sqr(point, node->val);
				if(dist_sqr < distance_sqaured)
					cluster.emplace_back(node);

				findNearestNodeClusterHelper(node->rightNode, point, distance, distance_sqaured, cluster, depth + 1);
				findNearestNodeClusterHelper(node->leftNode, point, distance, distance_sqaured, cluster, depth + 1);
				}
			*/

			/*
			if (test_if_the_hypersphere_and_hyperrectangle_intersect( point.coords, distance, node->coords_min, node->coords_max ))
				{

				if (test_if_the_hyperrectangle_is_fully_inside_the_sphere( point.coords, distance, node->coords_min, node->coords_max))
					{
					vector<std::shared_ptr<KdTreeNode>> subnodes = get_all_SubNodes(node);
					cluster.insert(cluster.end(), subnodes.begin(), subnodes.end());
					}

				else
					{
					float dist_sqr = calDist_sqr(point, node->val);
					if(dist_sqr < distance_sqaured)
						cluster.emplace_back(node);

					findNearestNodeClusterHelper(node->rightNode, point, distance, distance_sqaured, cluster, depth + 1);
					findNearestNodeClusterHelper(node->leftNode, point, distance, distance_sqaured, cluster, depth + 1);
					}
				}

			*/



			if (test_if_the_hypersphere_and_hyperrectangle_intersect( point.coords, distance, node->coords_min, node->coords_max ))
				{

				float dist_sqr = calDist_sqr(point, node->val);
				bool inside_sphere = false;
				if(dist_sqr < distance_sqaured)
					inside_sphere = true;

				if (inside_sphere && test_if_the_hyperrectangle_is_fully_inside_the_sphere( point.coords, distance, node->coords_min, node->coords_max))
					{
					vector <const KdTreeNode *> subnodes = get_all_SubNodes(node);
					cluster.insert(cluster.end(), subnodes.begin(), subnodes.end());
					}

				else
					{
					if (inside_sphere)
						cluster.emplace_back(node);

					findNearestNodeClusterHelper(node->rightNode, point, distance, distance_sqaured, cluster, depth + 1);
					findNearestNodeClusterHelper(node->leftNode, point, distance, distance_sqaured, cluster, depth + 1);
					}
				}


			/*
			if (test_if_the_hypersphere_and_hyperrectangle_intersect( point.coords, distance, node->coords_min, node->coords_max ))
				{
				float dist_sqr = calDist_sqr(point, node->val);
				if(dist_sqr < distance_sqaured)
					cluster.emplace_back(node);

				findNearestNodeClusterHelper(node->rightNode, point, distance, distance_sqaured, cluster, depth + 1);
				findNearestNodeClusterHelper(node->leftNode, point, distance, distance_sqaured, cluster, depth + 1);
				}

			*/


			/*
			int32_t curDim = depth % dimension;
			float dist_sqr = calDist_sqr(point, node->val);
			if(dist_sqr < distance_sqaured)
				cluster.emplace_back(node);

			float delta = point.coords[curDim] - node->val.coords[curDim];
			if(delta < 0)
				{
				findNearestNodeClusterHelper(node->leftNode, point, distance, distance_sqaured, cluster, depth + 1);
				if( -delta <= distance)
					//if (test_if_the_hypersphere_and_hyperrectangle_intersect( point.coords, distance, node->coords_min, node->coords_max ))
						{
						findNearestNodeClusterHelper(node->rightNode, point, distance, distance_sqaured, cluster, depth + 1);
						}
	 			}
			else
				{
				findNearestNodeClusterHelper(node->rightNode, point, distance, distance_sqaured, cluster, depth + 1);
				if(delta < distance)
					//if (test_if_the_hypersphere_and_hyperrectangle_intersect( point.coords, distance, node->coords_min, node->coords_max ))
						{
						findNearestNodeClusterHelper(node->leftNode, point, distance, distance_sqaured, cluster, depth + 1);
						}
				}


			*/
			}



		void write_KdTree_data_into_a_binary_data_stream(char *&data_stream_pointer, uint64_t &size) const
			{
			ERRORIF(data_stream_pointer != nullptr);

			uint64_t size_of_one_KdTreeNode = sizeof(PointType)*dimension * 3 + sizeof(size_t) + sizeof(uint8_t) + sizeof(uint8_t);

			cout << size_of_one_KdTreeNode << endl;

			size = sizeof(int32_t) + (uint64_t)count_number_of_nonnull_nodes()*size_of_one_KdTreeNode;

			data_stream_pointer = new char[size];

			// write dimension data
			uint64_t stream_position = 0;
			memcpy(&data_stream_pointer[stream_position],&dimension, sizeof(int32_t));
			stream_position+= sizeof(int32_t);

			write_KdTree_data_into_a_binary_data_stream_helper(root, data_stream_pointer, stream_position);

			ERRORIF(size != stream_position);

			}

		void write_KdTree_data_into_a_binary_data_stream_helper(const KdTreeNode * const node,  char * const &data_stream_pointer, uint64_t &stream_position) const
			{
			// write coordinate data
			memcpy(&data_stream_pointer[stream_position],&node->val.coords[0], sizeof(PointType)*dimension);
			stream_position+= sizeof(PointType)*dimension;
			// write index data
			memcpy(&data_stream_pointer[stream_position],&node->val.index, sizeof(size_t));
			stream_position+= sizeof(size_t);

			// write BB data
			memcpy(&data_stream_pointer[stream_position],&node->coords_max[0], sizeof(PointType)*dimension);
			stream_position+= sizeof(PointType)*dimension;
			memcpy(&data_stream_pointer[stream_position],&node->coords_min[0], sizeof(PointType)*dimension);
			stream_position+= sizeof(PointType)*dimension;

			// write sub-node data
			uint8_t rightnode_bool = 0;
			if(node->rightNode != nullptr)
				rightnode_bool = 1;
			uint8_t leftnode_bool = 0;
			if(node->leftNode != nullptr)
				leftnode_bool = 1;
			memcpy(&data_stream_pointer[stream_position],&rightnode_bool, sizeof(uint8_t));
			stream_position+= sizeof(uint8_t);
			memcpy(&data_stream_pointer[stream_position],&leftnode_bool, sizeof(uint8_t));
			stream_position+= sizeof(uint8_t);

			if(node->rightNode != nullptr)
				write_KdTree_data_into_a_binary_data_stream_helper(node->rightNode,  data_stream_pointer, stream_position);

			if(node->leftNode != nullptr)
				write_KdTree_data_into_a_binary_data_stream_helper(node->leftNode,  data_stream_pointer, stream_position);
			}


		void reconstruct_KdTree_data_from_a_binary_data_stream(const char * const &data_stream_pointer, const uint64_t size)
			{
			ERRORIF(data_stream_pointer == nullptr);
			ERRORIF(root != nullptr);

			// read dimension data
			uint64_t stream_position = 0;
			memcpy(&dimension, &data_stream_pointer[stream_position], sizeof(int32_t));
			stream_position+= sizeof(int32_t);

			reconstruct_KdTree_data_from_a_binary_data_stream_helper(root, data_stream_pointer, stream_position);

			ERRORIF(size != stream_position);

			// generate Bounding Box data
			/*vector <PointType> coords_min (dimension, 0);
			vector <PointType> coords_max (dimension, 0);
			add_Bounding_Box_data_to_the_tree_helper(root, coords_min, coords_max);*/
			}

		void reconstruct_KdTree_data_from_a_binary_data_stream_helper(KdTreeNode * & node,  const char * const &data_stream_pointer, uint64_t &stream_position)
			{
			node = new KdTreeNode();

			// read coordinate data
			node->val.coords = vector <PointType> (dimension,0);
			memcpy(&node->val.coords[0], &data_stream_pointer[stream_position], sizeof(PointType)*dimension);
			stream_position+= sizeof(PointType)*dimension;

			// read index data
			memcpy(&node->val.index, &data_stream_pointer[stream_position], sizeof(size_t));
			stream_position+= sizeof(size_t);

			// read BB data
			node->coords_max = vector <PointType> (dimension,0);
			memcpy(&node->coords_max[0], &data_stream_pointer[stream_position], sizeof(PointType)*dimension);
			stream_position+= sizeof(PointType)*dimension;
			node->coords_min = vector <PointType> (dimension,0);
			memcpy(&node->coords_min[0], &data_stream_pointer[stream_position], sizeof(PointType)*dimension);
			stream_position+= sizeof(PointType)*dimension;

			// read sub-node data
			node->rightNode = nullptr;
			node->leftNode = nullptr;

			uint8_t rightnode_bool;
			memcpy(&rightnode_bool, &data_stream_pointer[stream_position], sizeof(uint8_t));
			stream_position+= sizeof(uint8_t);
			uint8_t leftnode_bool;
			memcpy(&leftnode_bool, &data_stream_pointer[stream_position], sizeof(uint8_t));
			stream_position+= sizeof(uint8_t);

			if(rightnode_bool == 1)
				reconstruct_KdTree_data_from_a_binary_data_stream_helper(node->rightNode, data_stream_pointer, stream_position);

			if(leftnode_bool == 1)
				reconstruct_KdTree_data_from_a_binary_data_stream_helper(node->leftNode, data_stream_pointer, stream_position);
			}


		void write_KdTree_to_a_binary_file (const string fname) const
			{
			char *data_stream_pointer = nullptr;
			uint64_t size = 0;
			write_KdTree_data_into_a_binary_data_stream(data_stream_pointer, size);

			FILE * pFile;
			pFile = fopen ( fname.c_str() , "wb" );
			if (pFile==NULL)
				error(AT,FU, "Problems opening file: " + fname);

			size_t result = fwrite(data_stream_pointer,sizeof(char),size,pFile);
			if (result != size)
				error(AT,FU, "Problems writing kdtree to file: " + fname);
			fclose (pFile);

			delete [] data_stream_pointer;

			cout << "Writting kdtree to file: " << fname << endl;
			}

		void read_KdTree_from_a_binary_file (const string fname)
			{
			FILE * pFile;
			pFile = fopen ( fname.c_str() , "rb" );
			if (pFile==NULL)
				error(AT,FU, "Problems opening file: " + fname);

			fseek(pFile,0,SEEK_END);
			uint64_t size = ftell(pFile);
			rewind(pFile);

			char *data_stream_pointer = nullptr;
			data_stream_pointer = new char[size];

			size_t result = fread (data_stream_pointer,sizeof(char),size,pFile);
			if (result != size)
				error(AT,FU, "Problems reading kdtree from file: " + fname);
			fclose (pFile);

			reconstruct_KdTree_data_from_a_binary_data_stream(data_stream_pointer, size);

			delete [] data_stream_pointer;

			cout << "Reading kdtree from file: " << fname << endl;
			}

		};

}

