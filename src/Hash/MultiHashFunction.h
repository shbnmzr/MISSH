/*
 * MultiHashFunction.h
 *
 */

#ifndef HASH_MULTIHASHFUNCTION_H_
#define HASH_MULTIHASHFUNCTION_H_

#include "HashFunction.h"
#include <omp.h>
#include <math.h> 


inline static void GetHashes_speedup_multi_previous(const string& s_Str,
		const SpacedQmer_Multi& spaced_qmers,
		V_V_Hash_Err& vHashes,
		hash_type (*fConvertion)(char)) {

	auto get_hash = [&](size_t curr_spaced, size_t curr_idx_hash, const PreviusShift_Ext& curr_shift){
		Hash_Err& curr_hash = vHashes.get(curr_spaced, curr_idx_hash);

		if(curr_shift.current_sp_ptr->GetWeight() < curr_shift.GetSize())
			GetHash(s_Str, curr_idx_hash, *curr_shift.current_sp_ptr, curr_hash, fConvertion);
		else
		{
			size_t pos_hash_get = curr_idx_hash-curr_shift.shift_min;//la posizione dell'hash presa è la posizione attuale meno l'indice dello shift dove si fan meno cambiamenti
			const Hash_Err& prev_hash = vHashes.get(curr_shift.prev_qmer, pos_hash_get);

			compute_hash_for_speedup_previous(s_Str,
						curr_shift.current_sp_ptr->GetPosOne(), curr_shift.prev_sp_ptr->GetPosOne(),
						curr_shift,
						prev_hash,
						curr_idx_hash, curr_hash,
						fConvertion);
		}
	};

	//calcolo dimensione v_hash presupponendo che son tutti uguali
	long n_hashes = s_Str.size()-spaced_qmers[0].GetQ()+1;//Presumo lunghezza uguali tra i seed

	vHashes.clear();
	if(n_hashes>0)
	{
		const vector<PreviusShift_Ext_V>& v_shift_min = spaced_qmers.getShiftMin();
		vHashes.resize(spaced_qmers.size(), n_hashes, Hash_Err());

		for(size_t s = 0; s < spaced_qmers.size(); ++s)
		{
			if(!v_shift_min[s][0].isCorrectSpacedPrevious())
				GetHash(s_Str, 0, spaced_qmers[s], vHashes.get(s, 0), fConvertion);//primo da computare a parte
			else
				get_hash(s, 0, v_shift_min[s][0]);
		}

		//TODO: Attenzione a dimensione vettore hash e shift se di dimensioni diverse i seed
		////vHashes[0].size() == vHashes[1].size() == ... == vHashes[n].size()
		////v_shift_min[0].size() == v_shift_min[1].size() == ... == v_shift_min[n].size()
		size_t lim_max = vHashes.size_row(0);
		size_t lim_min = v_shift_min[0].size() < lim_max ? v_shift_min[0].size() : lim_max;

		for(size_t i = 1; i < lim_min; ++i)
			for(size_t ss = 0; ss < spaced_qmers.size(); ++ss)
				get_hash(ss, i, v_shift_min[ss][i]);
		for(size_t i = lim_min; i < lim_max; ++i)//vHashes[0].size() == vHashes[1].size() == ... == vHashes[n].size()
			for(size_t ss = 0; ss < spaced_qmers.size(); ++ss)
				get_hash(ss, i, v_shift_min[ss].back());
	}
	else
		vHashes.resize(spaced_qmers.size(), 0, Hash_Err());
}

inline static void GetHashes_speedup_multi_previous_Rotated(const string& s_Str,
		const SpacedQmer_Multi& spaced_qmers,
		V_V_Hash_Err& vHashes,
		hash_type (*fConvertion)(char)) {

	auto get_hash = [&](size_t curr_spaced, size_t curr_idx_hash, const PreviusShift_Ext& curr_shift){
		Hash_Err& curr_hash = vHashes.get(curr_spaced, curr_idx_hash);

		if(curr_shift.current_sp_ptr->GetWeight() < curr_shift.GetSize())
			GetHash(s_Str, curr_idx_hash, *curr_shift.current_sp_ptr, curr_hash, fConvertion);
		else
		{
			size_t pos_hash_get = curr_idx_hash-curr_shift.shift_min;//la posizione dell'hash presa è la posizione attuale meno l'indice dello shift dove si fan meno cambiamenti
			const Hash_Err& prev_hash = vHashes.get(curr_shift.prev_qmer, pos_hash_get);

			compute_hash_for_speedup_previous(s_Str,
						curr_shift.current_sp_ptr->GetPosOne(), curr_shift.prev_sp_ptr->GetPosOne(),
						curr_shift,
						prev_hash,
						curr_idx_hash, curr_hash,
						fConvertion);
		}
	};

	//calcolo dimensione v_hash presupponendo che son tutti uguali
	long n_hashes = s_Str.size()-spaced_qmers[0].GetQ()+1;//Presumo lunghezza uguali tra i seed

	vHashes.clear();
	if(n_hashes>0)
	{
		const vector<PreviusShift_Ext_V>& v_shift_min = spaced_qmers.getShiftMinRotated();
		vHashes.resize(spaced_qmers.size(), n_hashes, Hash_Err());

		const PreviusShift_Ext_V& zero_prev = v_shift_min[0];
		for(size_t s = 0; s < spaced_qmers.size(); ++s)
		{
			const PreviusShift_Ext& zero = zero_prev[s];
			if(!zero.isCorrectSpacedPrevious())
				GetHash(s_Str, 0, *zero.current_sp_ptr, vHashes.get(s, 0), fConvertion);//primo da computare a parte
			else
				get_hash(s, 0, zero);
		}

		//TODO: Attenzione a dimensione vettore hash e shift se di dimensioni diverse i seed
		////vHashes[0].size() == vHashes[1].size() == ... == vHashes[n].size()
		////v_shift_min[0].size() == v_shift_min[1].size() == ... == v_shift_min[n].size()
		size_t lim_max = vHashes.size_row(0);
		size_t lim_min = v_shift_min.size() < lim_max ? v_shift_min.size() : lim_max;

		for(size_t i = 1; i < lim_min; ++i)
		{
			const PreviusShift_Ext_V& prev_i = v_shift_min[i];
			for(size_t ss = 0; ss < spaced_qmers.size(); ++ss)
				get_hash(ss, i, prev_i[ss]);
		}
		for(size_t i = lim_min; i < lim_max; ++i)//vHashes[0].size() == vHashes[1].size() == ... == vHashes[n].size()
		{
			const PreviusShift_Ext_V& prev_i = v_shift_min.back();
			for(size_t ss = 0; ss < spaced_qmers.size(); ++ss)
				get_hash(ss, i, prev_i[ss]);
		}
	}
	else
		vHashes.resize(spaced_qmers.size(), 0, Hash_Err());
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// ISSH multiseed


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/**
 * Compute the hashing for a column. All the previous position are computed based on the ISSH single information (for an hash relative to a single seed only parts of hashes of the same seed are used), 
 * but the computation of the new character is done only once - using the encode function - and than used for all the hashes in the column.
 * 
 * @param curr_idx_hash index of the column of hashes that will be computed
 * @param VV_shifts it contains the data structures needed for the ISSH computation, for both standard and transient parts
 * @param vvHash structure that will contain all the computed hashes, each row correspond to a different spaced seed
 * @param v_pos_one array containing all the positions in which the ones are present inside each spaced seed
 * @param s_Str string containing the sequence of DNA
 * @param fConvertion encode function
*/
inline static void compute_hash_with_ISSH_multi_v1(			size_t curr_idx_hash,
															const vector<V_V_PreviusShift>& VV_shifts,
															Hash_Err_V_V& vvHash,
															const vector<Position>& v_pos_one,
															const string& s_Str,
															hash_type (*fConvertion)(char))
{

	// At this point the transitory part is done, in each position it is possible to retrive all the 
	// positions, but the last one.
	// We know the last position and therefore the corresponding value inside the DNA sequence.
	const Position& pos_not_covered_yet = VV_shifts[0][0][0].one_to_change;
	const size_t& i_to_change = pos_not_covered_yet[0];
	
	// identify the position of the unknown character inside the DNA sequence
	size_t index_char = curr_idx_hash+v_pos_one[0][i_to_change];
	hash_type ch = (*fConvertion)(s_Str[index_char]);
	
	for(size_t j = 0; j < VV_shifts.size(); ++j)
	{
		Hash_Err& curr_hash = vvHash[j][curr_idx_hash];
		curr_hash.hash = 0;
		
		for(size_t k = 0; k < VV_shifts[j][0].size(); k++)
		{
			// get the previous has from which to recover some positions 
			size_t pos_hash_get = curr_idx_hash-VV_shifts[j][0][k].shift_min;
			const Hash_Err& prev_hash = vvHash[j][pos_hash_get];
			
			hash_type partial_hash = prev_hash.hash;
			const PreviousShift& curr_sp_shift= VV_shifts[j][0][k];
			// shift the hash by the precomputed amount
			partial_hash >>= 2*curr_sp_shift.one_exit;

			// clear all the unwanted information
			partial_hash &= curr_sp_shift.mask;

			// manage the errors TODO: probably it needs fixing
			if(!prev_hash.isCorrect())
			{
				long curr_pos_one = 0;
				for(size_t e = 0; e < prev_hash.size_error(); ++e)
					if((curr_pos_one = prev_hash[e]-curr_sp_shift.one_exit) >= 0)
						if(v_pos_one[j][prev_hash[e]]-curr_sp_shift.shift_min == v_pos_one[j][curr_pos_one])
							curr_hash.push_back_error(curr_pos_one);
			}
			curr_hash.hash |= partial_hash; //add together the information to compute the hash

		}
		if(ch == 4) //convertion error
			curr_hash.push_back_error(i_to_change);
		else
			curr_hash.hash |= ch << (i_to_change * 2); //add the code for the new character

	}
};




/**
 * Compute the hashing for the given sequence based on the group of spaced seeds. It computes the transitory part for each seeds separately
 * and than starts computing a column at a time using the same information as ISSH, but using the encode function for the new character only 
 * once for all the hashes in a column
 * 
 * @param s_Str string containing the sequence of DNA
 * @param spaced_qmers array containing the information about each seed
 * @param VV_shifts it contains the data structures needed for the ISSH computation, for both standard and transient parts
 * @param v_pos_one array containing all the positions in which the ones are present inside each spaced seed
 * @param max_transient_length maximum overall length of the transient 
 * @param vvHash structure that will contain all the computed hashes, each row correspond to a different spaced seed
 * @param fConvertion encode function
*/
inline static void GetHashes_with_ISSH_multi_v1(    const string& s_Str,
										 			const vector<SpacedQmer>& spaced_qmers,
													const vector<V_V_PreviusShift>& VV_shifts,
													const vector<Position>& v_pos_one,
													const size_t max_transient_length,
													Hash_Err_V_V& vvHash,
										 			hash_type (*fConvertion)(char)
												)
{
	// Compute the number of hashes that needs to be computed
	long n_hashes = s_Str.size()-spaced_qmers[0].GetQ()+1;
	
	if(n_hashes>0)
	{
		size_t i=1;
		size_t shift_group;
		size_t curr_max;
		// for each of the spaced seeds initialize result struct and compute transient part
		for(size_t j = 0; j < spaced_qmers.size(); ++j)
		{	
			// clear the struct that will contain the hashing results
			vvHash[j].clear();
			vvHash[j].resize(n_hashes);

			const V_V_PreviusShift& curr_seed = VV_shifts[j];
			curr_max = curr_seed.size();
			// the first hash needs to be computed anew for each seed
			GetHash(s_Str, 0, spaced_qmers[j], vvHash[j][0], fConvertion);
			i=1;
			
			for(; i<(size_t)n_hashes && i<max_transient_length; i++)
			{
				// if there's no recepy for what previous hash to use compute it anew 
				if(curr_seed[i].size()==0)
				{
					GetHash(s_Str, i, spaced_qmers[j], vvHash[j][i], fConvertion);
				}
				// otherwise use ISSH computation
				else
				{	
					shift_group = i >= curr_max ? 0 : i;
					compute_hash_with_ISSH(i, curr_seed[shift_group], vvHash[j], v_pos_one[j], s_Str, fConvertion);
				}
			}
		}
		// once the transient part is computed the standard computation can start for all the remaining columns
		for(; i < (size_t)n_hashes; i++)
		{
			//compute the hashes at position i for all the seeds 
			compute_hash_with_ISSH_multi_v1(i, VV_shifts, vvHash, v_pos_one, s_Str, fConvertion);
		
		}
	}
	// return an empty struct if no hash can be computed
	else
	{
		// clear hashing struct
		for(size_t j = 0; j < VV_shifts.size(); j++)
		{
			vvHash[j].clear();
		}
	}
}



/**
 * Compute the hashing for a column. For each hash the information about the hashes in the previous columns and in the same
 * column for the rows above can be used.
 * 
 * @param curr_idx_hash index of the column of hashes that will be computed
 * @param VV_shifts it contains the data structures needed for the ISSH computation, for both standard and transient parts
 * @param vvHash structure that will contain all the computed hashes, each row correspond to a different spaced seed
 * @param s_Str string containing the sequence of DNA
 * @param fConvertion encode function
*/
inline static void compute_hash_with_ISSH_multi_col(		size_t curr_idx_hash,
															const MultiSeedInfo& VV_shifts,
															Hash_Err_V_V& vvHash,
															const string& s_Str,
															hash_type (*fConvertion)(char))
{

	// find what is the current situation, if the computation is standard or it is a transient step
	// the seed length is considered by 1+last position in which a one is present in the seed
	// since n_hashes = sequence_length - seed_length + 1 we can compute it as
	
	int shift_group = curr_idx_hash < VV_shifts[0].group_previous.size()-1 ? curr_idx_hash + 1 : 0;
	// for each of the seed in the group	
	for(size_t j = 0; j < VV_shifts.size(); ++j)
	{
		// Get the reference to the current hash to compute
		Hash_Err& curr_hash = vvHash[j][curr_idx_hash];
		// Initialize it to zero
		curr_hash.hash = 0;
		const groupPrevious& curr_group = VV_shifts[j].group_previous[shift_group];
		const V_PreviusShiftMulti& curr_group_prev = curr_group.prev;
		
		// for faster performance the first hash is computed by calling GetHash instead of recovering all the positions later
		if(VV_shifts[j].group_previous[shift_group].prev.size()==0)
		{	
			GetHashFromPosOne(s_Str, curr_idx_hash, VV_shifts[j].pos_ones, vvHash[j][curr_idx_hash], fConvertion);
		}
		else
		{	
			for(size_t k = 0; k < curr_group_prev.size(); k++)
			{
				const PreviousShiftMulti& curr_sp_shift= curr_group_prev[k];
				// identify position of the previous hash
				size_t pos_hash_get = curr_idx_hash-curr_sp_shift.shift_min;
				// get the previous hash
				const Hash_Err& prev_hash = vvHash[curr_sp_shift.seed_num][pos_hash_get];
				// take a copy of the previous hash
				hash_type partial_hash = prev_hash.hash;
				// shift the hash based on the offset
				if(curr_sp_shift.one_exit>=0)
					partial_hash >>= 2*curr_sp_shift.one_exit;
				else
					partial_hash <<= -2*curr_sp_shift.one_exit;
				// clear useless position by using the mask
				partial_hash &= curr_sp_shift.mask;
				
				// Error positions in the previous hash should be considered as error even in the new hash 
				if(!prev_hash.isCorrect())
				{
					//TODO: manage errors in the sequence
				}
				curr_hash.hash |= partial_hash;
			}

			// compute from the seqence all the missing positions
			const Position& pos_not_covered_yet = curr_group.not_covered;
			const Position& pos_one = VV_shifts[j].pos_ones;
			for(size_t i = 0; i < pos_not_covered_yet.size(); i++)
			{
				size_t i_to_change = pos_not_covered_yet[i];
				hash_type ch = (*fConvertion)(s_Str[curr_idx_hash+pos_one[i_to_change]]);
				if(ch == 4) //Errore conversione
					curr_hash.push_back_error(i_to_change);
				else
					curr_hash.hash |= ch << (i_to_change * 2);//OR, add right values on right positions
			}
		}
	}
};





/**
 * Compute the hashing for the given sequence based on the group of spaced seeds. It computes a column at a time by using all the 
 * possible previously computed hashes for recovering positions
 * 
 * @param s_Str string containing the sequence of DNA
 * @param VV_shifts it contains the data structures needed for the ISSH computation by column
 * @param vvHash structure that will contain all the computed hashes, each row correspond to a different spaced seed
 * @param fConvertion encode function
*/
inline static void GetHashes_with_ISSH_multi_col(   const string& s_Str,
													const MultiSeedInfo& VV_shifts,
													Hash_Err_V_V& vvHash,
										 			hash_type (*fConvertion)(char)
	 											)
{
	// Total number of hashes that can be computed
	long n_hashes = s_Str.size()-VV_shifts[0].pos_ones[VV_shifts[0].pos_ones.size()-1];
	if(n_hashes>0)
	{
		// initialize the matrix of hashes, one row for each seed
		for(size_t j = 0; j < VV_shifts.size(); j++)
		{	
			vvHash[j].clear();
			vvHash[j].resize(n_hashes);
		}
		// compute a column at a time
		for(size_t i = 0; i < (size_t)n_hashes; i++)
		{
			compute_hash_with_ISSH_multi_col(i, VV_shifts, vvHash, s_Str, fConvertion);
		}
	}else
	{
		// clear hashing struct
		for(size_t j = 0; j < VV_shifts.size(); j++)
		{
			vvHash[j].clear();
		}
	}
}

/**
 * Compute the hashing for a column. For each hash the information about the hashes in the previous columns and in the same
 * column for the rows above can be used.
 * 
 * @param curr_idx_hash index of the column of hashes that will be computed
 * @param offset it tells which step of the transient we are considering since multiple transients are needed
 * @param VV_shifts it contains the data structures needed for the ISSH computation, for both standard and transient parts
 * @param vvHash structure that will contain all the computed hashes, each row correspond to a different spaced seed
 * @param s_Str string containing the sequence of DNA
 * @param fConvertion encode function
*/
inline static void compute_hash_with_ISSH_multi_col_parallel(size_t curr_idx_hash,
															size_t offset,
															const MultiSeedInfo& VV_shifts,
															Hash_Err_V_V& vvHash,
															const string& s_Str,
															hash_type (*fConvertion)(char))
{
	

	// find what is the current situation, if the computation is standard or it is a transient step
	int shift_group = offset < VV_shifts[0].group_previous.size()-1 ? offset + 1 : 0;
	
	for(size_t j = 0; j < VV_shifts.size(); ++j)
	{
		// Get the reference to the current hash to compute
		
		Hash_Err& curr_hash = vvHash[j][curr_idx_hash];
	
		// Initialize it to zero
		curr_hash.hash = 0;
		const groupPrevious& curr_group = VV_shifts[j].group_previous[shift_group];
		const V_PreviusShiftMulti& curr_group_prev = curr_group.prev;
		
		// for faster performance the first hash is computed by calling GetHash instead of recovering all the positions later
		if(VV_shifts[j].group_previous[shift_group].prev.size()==0)
		{	
			GetHashFromPosOne(s_Str, curr_idx_hash, VV_shifts[j].pos_ones, vvHash[j][curr_idx_hash], fConvertion);
		}
		else
		{	
			for(size_t k = 0; k < curr_group_prev.size(); k++)
			{
				const PreviousShiftMulti& curr_sp_shift= curr_group_prev[k];
				// identify position of the previous hash
				size_t pos_hash_get = curr_idx_hash-curr_sp_shift.shift_min;
				// get the previous hash
				const Hash_Err& prev_hash = vvHash[curr_sp_shift.seed_num][pos_hash_get];
				
				// take a copy of the previous hash
				hash_type partial_hash = prev_hash.hash;
				// shift the hash based on the offset
				if(curr_sp_shift.one_exit >= 0)
					partial_hash >>= 2*curr_sp_shift.one_exit;
				else
					partial_hash <<= -2*curr_sp_shift.one_exit;
				// clear useless position by using the mask
				partial_hash &= curr_sp_shift.mask;
				
				// Error positions in the previous hash should be considered as error even in the new hash 
				if(!prev_hash.isCorrect())
				{
					//TODO: manage errors in the sequence
				}
				curr_hash.hash |= partial_hash;
			}

			// compute from the seqence all the missing positions
			const Position& pos_not_covered_yet = curr_group.not_covered;
			const Position& pos_one = VV_shifts[j].pos_ones;
			for(size_t i = 0; i < pos_not_covered_yet.size(); i++)
			{
				size_t i_to_change = pos_not_covered_yet[i];
				hash_type ch = (*fConvertion)(s_Str[curr_idx_hash+pos_one[i_to_change]]);
				if(ch == 4) //Errore conversione
					curr_hash.push_back_error(i_to_change);
				else
					curr_hash.hash |= ch << (i_to_change * 2);//OR, add right values on right positions
			}
		}
	}
};



/**
 * Compute the hashing for the given sequence based on the group of spaced seeds. It computes a column at a time by using all the 
 * possible previously computed hashes for recovering positions. The hashing is divided into blocks that can be computed by different threads.
 * Each thread needs to start the computation with a transient part.
 * 
 * @param s_Str string containing the sequence of DNA
 * @param VV_shifts it contains the data structures needed for the ISSH computation by column
 * @param vvHash structure that will contain all the computed hashes, each row correspond to a different spaced seed
 * @param fConvertion encode function
*/
inline static void GetHashes_with_ISSH_multi_col_parallel
												(   const string& s_Str,
													const MultiSeedInfo& VV_shifts,
													Hash_Err_V_V& vvHash,
										 			hash_type (*fConvertion)(char)
	 											)
{
	// Total number of hashes that we want to compute
	// the seed length is considered by 1+last position in which a one is present in the seed
	// since n_hashes = sequence_length - seed_length + 1 we can compute it as
	
	long n_hashes = s_Str.size()-VV_shifts[0].pos_ones[VV_shifts[0].pos_ones.size()-1];
	// clear hashing struct
	
	if(n_hashes>0)
	{
		// initialize the matrix of hashes, one for each seed
		for(size_t j = 0; j < VV_shifts.size(); j++)
		{	
			vvHash[j].clear();
			vvHash[j].resize(n_hashes);
		}
		// compute a column at a time
		
		int offset; // declared as a private variable for each thread 
		#pragma omp parallel private(offset)// num_threads(thread_num) 
		{
			int thread_num = omp_get_num_threads();
			offset = 0; // declared as a private variable for each thread

			#pragma omp for schedule(static, (int)ceil((n_hashes+0.0)/thread_num))
			for(size_t i = 0; i < (size_t)n_hashes; i++)
			{	
				compute_hash_with_ISSH_multi_col_parallel(i, offset, VV_shifts, vvHash, s_Str, fConvertion);
				offset++;
			}
		}
	}
	else
	{
		for(size_t j = 0; j < VV_shifts.size(); j++)
		{
			vvHash[j].clear();
		}
	}
}






/**
 * Compute the hashing for the given sequence based on the group of spaced seeds. It computes a row at a time by using all the 
 * possible previously computed hashes for recovering positions (previous rows and same row but previous column). 
 * 
 * @param s_Str string containing the sequence of DNA
 * @param rowInfo it contains the data structures needed for the ISSH computation by row
 * @param vvHash structure that will contain all the computed hashes, each row correspond to a different spaced seed
 * @param fConvertion encode function
*/
inline static void GetHashes_with_ISSH_multi_row(   const string& s_Str,
													const MultiSeedInfoRow& rowInfo,
													Hash_Err_V_V& vvHash,
										 			hash_type (*fConvertion)(char)
	 											)
{
	const MultiSeedInfo & VV_shifts = rowInfo.info;	
	// Total number of hashes that we want to compute
	// the seed length is considered by 1+last position in which a one is present in the seed
	// since n_hashes = sequence_length - seed_length + 1 we can compute it as
	long n_hashes = s_Str.size()-VV_shifts[0].pos_ones[VV_shifts[0].pos_ones.size()-1];
	
	size_t transient1 = rowInfo.transient1;
	size_t transient2 = rowInfo.transient2;
	
	if(n_hashes>0 && n_hashes>=(long)(transient1+transient2))
	{
		// initialize the matrix of hashes, one for each seed
		for(size_t j = 0; j < VV_shifts.size(); j++)
		{
			// clear hashing struct
			vvHash[j].clear();
			vvHash[j].resize(n_hashes);
			
			const V_groupPrevious & V_shift = VV_shifts[j].group_previous;
			// compute all the hashes for the current spaced seed
			for(size_t curr_idx_hash = 0; curr_idx_hash < (size_t)n_hashes; curr_idx_hash++)
			{
				int shift_group = (curr_idx_hash >= transient1 && curr_idx_hash < n_hashes - transient2) ? 0 : ((curr_idx_hash < transient1) ? (curr_idx_hash + 1) :(V_shift.size() - (n_hashes - curr_idx_hash)));
				const groupPrevious& curr_group = V_shift[shift_group];
				const V_PreviusShiftMulti& curr_group_prev = curr_group.prev;
			
				// Get the reference to the current hash to compute
				Hash_Err& curr_hash = vvHash[j][curr_idx_hash];
				// Initialize it to zero
				curr_hash.hash = 0;
				//cout<<"Curr hash = "<<vvHash[j][curr_idx_hash].hash<<endl;
				// for faster performance the first hash is computed by calling GetHash instead of recovering all the positions later
				if(VV_shifts[j].group_previous[shift_group].prev.size()==0)
				{	
					GetHashFromPosOne(s_Str, curr_idx_hash, VV_shifts[j].pos_ones, vvHash[j][curr_idx_hash], fConvertion);
					//cout << vvHash[j][curr_idx_hash].hash<<endl;
				}
				else
				{	
					// use previously computed hashes to recover the position
					for(size_t k = 0; k < curr_group_prev.size(); k++)
					{	
						const PreviousShiftMulti& curr_sp_shift= curr_group_prev[k];
						// identify position of the previous hash
						size_t pos_hash_get = curr_idx_hash-curr_sp_shift.shift_min;
						// get the previous hash
						const Hash_Err& prev_hash = vvHash[curr_sp_shift.seed_num][pos_hash_get];
						// take a copy of the previous hash
						hash_type partial_hash = prev_hash.hash;
						// shift the hash based on the offset
						if(curr_sp_shift.one_exit>=0)
							partial_hash >>= 2*curr_sp_shift.one_exit;
						else
							partial_hash <<= -2*curr_sp_shift.one_exit;
						// clear useless position by using the mask
						partial_hash &= curr_sp_shift.mask;
						
						// Error positions in the previous hash should be considered as error even in the new hash 
						if(!prev_hash.isCorrect())
						{
							//TODO: manage errors in the sequence
						}
						curr_hash.hash |= partial_hash;
					}

					// compute from the sequence all the missing positions
					const Position& pos_not_covered_yet = curr_group.not_covered;
					const Position& pos_one = VV_shifts[j].pos_ones;
					for(size_t i = 0; i < pos_not_covered_yet.size(); i++)
					{
						size_t i_to_change = pos_not_covered_yet[i];
						hash_type ch = (*fConvertion)(s_Str[curr_idx_hash+pos_one[i_to_change]]);
						if(ch == 4) //Errore conversione
							curr_hash.push_back_error(i_to_change);
						else
							curr_hash.hash |= ch << (i_to_change * 2);//OR, add right values on right positions
					}
				}
			}
		}
	}else if(n_hashes<(long)(transient1+transient2) && n_hashes>0)
	{	
		for(size_t j = 0; j < VV_shifts.size(); j++)
		{
			// clear hashing struct
			vvHash[j].clear();
			vvHash[j].resize(n_hashes);
			for(size_t curr_idx_hash = 0; curr_idx_hash < (size_t)n_hashes; curr_idx_hash++)
			{
				GetHashFromPosOne(s_Str, curr_idx_hash, VV_shifts[j].pos_ones, vvHash[j][curr_idx_hash], fConvertion);
			}
		}	
	}else
	{
		for(size_t j = 0; j < VV_shifts.size(); j++)
		{
			vvHash[j].clear();
		}
	}
}





#endif /* HASH_MULTIHASHFUNCTION_H_ */
