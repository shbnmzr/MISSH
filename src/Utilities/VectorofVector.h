/*
 * VectorofVector.h
 *
 *  Created on: 21/feb/2017
 *      Author: samuele
 */

#ifndef UTILITIES_VECTOROFVECTOR_H_
#define UTILITIES_VECTOROFVECTOR_H_
#include <vector>
#include <memory>
#include <stdexcept>
#include <cassert>
#include <numeric>
using namespace std;
#define NDEBUG

template<typename T>
class Vector_of_Vector {
	struct Row {
		size_t start;
		size_t size;
	};
public:
	typedef T value_type;
	typedef reference_wrapper<value_type> reference_value_type;

	Vector_of_Vector() {

	}
	//init matrix
	Vector_of_Vector(size_t rows, size_t cols) {
		this->resize(rows, cols);
	}
	Vector_of_Vector(size_t rows, size_t cols, const value_type& val) {
		this->resize(rows, cols, val);
	}
	//init vector of vector with different size
	Vector_of_Vector(const vector<size_t> size_rows){
		this->resize(size_rows);
	}
	Vector_of_Vector(const vector<size_t> size_rows, const value_type& val) {
		this->resize(size_rows, val);
	}
	virtual ~Vector_of_Vector() {

	}
	inline const value_type& get(size_t row, size_t col) const{
//		if(col >= this->size_row(row))
//			throw out_of_range("out_of_range error: col >= this->size_row(row) -> " + to_string(col) + " >= " + to_string(this->size_row(row)));
		return this->data[this->rows[row].start+col];
	}
	inline value_type& get(size_t row, size_t col) {
//		if(col >= this->size_row(row))
//			throw out_of_range("out_of_range error: col >= this->size_row(row) -> " + to_string(col) + " >= " + to_string(this->size_row(row)));
		return this->data[this->rows[row].start+col];
	}
	//return only reference for object in matrix,
	//insert in this vector not insert in matrix,
	//for insert use insert method
	inline void get(size_t row, vector<reference_value_type>& result) {
		if(row >= this->rows.size())
			throw out_of_range("out_of_range error: row >= this->rows.size() -> " +  to_string(row) + " >= " +  to_string(this->rows.size()));
		size_t row_size = this->size_row(row);
		result.clear();
		result.reserve(row_size);
		for(size_t i = this->rows[row].start; i < this->rows[row].start + row_size; ++i)
			result.push_back(this->data[i]);
	}
	inline void insert(size_t row, size_t col, const T& value) {
		if(row > this->rows.size())
			throw out_of_range("out_of_range error: row > this->rows.size() -> " + to_string(row) + " > " + to_string(this->rows.size()));
		if(col > (row == this->rows.size() ? 0:size_row(row)))
			throw out_of_range("out_of_range error: col > (row == this->rows.size() ? 0:size_row(row)) -> " + to_string(col) + " > " + to_string((row == this->rows.size() ? 0:size_row(row))));
		auto pos_insert = advance(this->data.begin(), row*col);
		this->data.insert(pos_insert, value);
		if(row == this->rows.size())
		{
			Row new_row;
			new_row.start = this->data.size()-1;
			new_row.size = 1;
			this->rows.push_back(move(new_row));
		}
		else
		{
			++this->rows[row].size;
			//#pragma omp parallel for
			for(size_t i = row+1; i < this->rows.size(); ++i)
				++this->rows[i].start;
		}
	}
	inline void insert(size_t row, size_t col, T&& value) {
		if(row > this->rows.size())
			throw out_of_range("out_of_range error: row > this->rows.size() -> " + to_string(row) + " > " + to_string(this->rows.size()));
		if(col > (row == this->rows.size() ? 0:size_row(row)))
			throw out_of_range("out_of_range error: col > (row == this->rows.size() ? 0:size_row(row)) -> " + to_string(col) + " > " + to_string((row == this->rows.size() ? 0:size_row(row))));
		auto pos_insert = advance(this->data.begin(), row*col);
		this->data.insert(pos_insert, value);
		if(row == this->rows.size())
		{
			Row new_row;
			new_row.start = this->data.size()-1;
			new_row.size = 1;
			this->rows.push_back(move(new_row));
		}
		else
		{
			++this->rows[row].size;
			//#pragma omp parallel for
			for(size_t i = row+1; i < this->rows.size(); ++i)
				++this->rows[i].start;
		}
	}
	inline size_t size_row(size_t row) const {
//		if(row >= this->rows.size())
//			throw out_of_range("out_of_range error: row >= this->rows.size() -> " +  to_string(row) + " >= " +  to_string(this->rows.size()));
		return this->rows[row].size;
	}
	inline size_t size() const {
		return this->data.size();
	}
	inline void clear() {
		this->data.clear();
		this->rows.clear();
	}
	inline void swap(Vector_of_Vector& swap_el) {
		this->data.swap(swap_el.data);
		this->rows.swap(swap_el.rows);
	}
	inline void shrink_to_fit() {
		this->data.shrink_to_fit();
		this->rows.shrink_to_fit();
	}
	inline void reserve(size_t rows, size_t cols) {
		this->data.reserve(rows*cols);
		this->rows.reserve(rows);
	}
	inline void resize(size_t rows, size_t cols) {
		this->resize(rows, cols, value_type());
	}
	inline void resize(size_t rows, size_t cols, const value_type& val) {
		this->data.resize(rows*cols, val);
		this->rows.resize(rows);
		//#pragma omp parallel for
		for(size_t i = 0; i < this->rows.size(); ++i)
		{
			this->rows[i].start = i*cols;
			this->rows[i].size = cols;
		}
	}
	inline void resize(const vector<size_t> size_rows){
		this->resize(size_rows, value_type());
	}
	inline void resize(const vector<size_t> size_rows, const value_type& val) {
		size_t tot = accumulate(size_rows.begin(), size_rows.end(), 0);
		this->data.resize(tot, val);
		this->rows.resize(size_rows.size());
		//#pragma omp parallel for
		for(size_t i = 0; i < this->rows.size(); ++i)
			this->rows[i].size = size_rows[i];
		size_t start_idx = 0;
		for(size_t i = 0; i < this->rows.size(); ++i)
		{
			this->rows[i].start = start_idx;
			start_idx += size_rows[i];
		}
	}

private:
	vector<value_type> data;
	vector<Row> rows;
};

#endif /* UTILITIES_VECTOROFVECTOR_H_ */
