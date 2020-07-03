/*
 * HashType.h
 *
 *  Created on: 10/feb/2017
 *      Author: samuele
 */

#ifndef HASH_HASHTYPE_H_
#define HASH_HASHTYPE_H_

#include <vector>
#include <cstdint>
#include <memory>
#include "../Spaced/SpacedQmer.h"
#include "../Utilities/VectorofVector.h"
using namespace std;

//#include <boost/multiprecision/cpp_int.hpp>
//typedef uint128_t hash_type;
typedef uint64_t hash_type;

struct Hash_Err {
public:
	typedef shared_ptr<Position> Position_Ptr;
private:
	Position_Ptr err_pos;//indice dell'1 dello spaced seed che contiene l'errore
public:
	hash_type hash;

	Hash_Err(){
		this->hash = 0;
	}
	inline void reset() {
		this->hash = 0;
		this->err_pos.reset();
	}

	inline void sub_pos_err(size_t num_pos_to_subtract) {
		if(this->err_pos)
		{
			Position_Ptr tmp_err = make_shared<Position>();
			tmp_err->reserve(this->err_pos->size()*2);
			long curr_pos_one = 0;
			for(const size_t& err : *this->err_pos)
				if((curr_pos_one = err-num_pos_to_subtract) >= 0)
					tmp_err->push_back(curr_pos_one);
			this->err_pos = tmp_err;
		}
	}
	inline void add_pos_err(size_t num_pos_to_add) {
		if(this->err_pos)
		{
			Position_Ptr tmp_err = make_shared<Position>();
			tmp_err->reserve(this->err_pos->size()*2);
			for(const size_t& err : *this->err_pos)
				tmp_err->push_back(err+num_pos_to_add);
			this->err_pos = tmp_err;
		}
	}
	inline void sub_pos_err(size_t num_pos_to_subtract, const Hash_Err& hash_err) {
		if(hash_err.err_pos)
		{
			if(!this->err_pos)
				this->err_pos = make_shared<Position>();

			long curr_pos_one = 0;
			for(const size_t& err : *hash_err.err_pos)
				if((curr_pos_one = err-num_pos_to_subtract) >= 0)
					this->err_pos->push_back(curr_pos_one);
		}
	}
	inline void add_pos_err(size_t num_pos_to_add, const Hash_Err& hash_err) {
		if(hash_err.err_pos)
		{
			if(!this->err_pos)
				this->err_pos = make_shared<Position>();
			for(const size_t& err : *hash_err.err_pos)
				this->err_pos->push_back(err+num_pos_to_add);
		}
	}
	inline void sort_uniq_err() {
		if(this->err_pos)
		{
			sort(this->err_pos->begin(), this->err_pos->end());
			auto it = unique(this->err_pos->begin(), this->err_pos->end());
			this->err_pos->resize(distance(this->err_pos->begin(), it));
		}
	}

	inline bool isCorrect() const{
		if(this->err_pos)
			return this->err_pos->empty();
		else
			return !this->err_pos;
	};
	inline void create_error() {
		if(!this->err_pos)
			this->err_pos = make_shared<Position>();
	}
	inline void push_back_error(const Position::value_type& val) {
		if(!this->err_pos)
			this->err_pos = make_shared<Position>();
		this->err_pos->push_back(val);
	}
	inline size_t size_error() const {
		if(this->err_pos)
				return this->err_pos->size();
			else
				return 0;
	}
	inline const Position::value_type& operator[] (size_t i) const {return (*this->err_pos)[i];}
	inline Position::value_type& operator[] (size_t i) {return (*this->err_pos)[i];}
};
typedef vector<Hash_Err> Hash_Err_V;
typedef vector<Hash_Err_V> Hash_Err_V_V;
typedef Vector_of_Vector<Hash_Err> V_V_Hash_Err;

#endif /* HASH_HASHTYPE_H_ */
