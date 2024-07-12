/*
 * HashFunction.h
 *
 */

#ifndef HASHFUNCTION_H_
#define HASHFUNCTION_H_

#include "HashType.h"
#include "../Spaced/SpacedQmer_Multi.h"
#include "../Spaced/MultiSpacedQmer.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <bitset>

inline static hash_type CharToInt(char ch)
{
	if(ch == 'A')
		return 0;
	if(ch == 'C')
		return 1;
	if(ch == 'G')
		return 2;
	if(ch == 'T')
		return 3;
	return 4; //ERROR CODE
}

inline static void printhashes(Hash_Err_V& v_hash_err)
{
	for (size_t i=0; i< v_hash_err.size(); i++)
	cout<< v_hash_err[i].hash;
}

inline static hash_type CharToIntComplement(char ch)
{
	if(ch == 'A')
		return 3;
	if(ch == 'C')
		return 2;
	if(ch == 'G')
		return 1;
	if(ch == 'T')
		return 0;
	return 4; //ERROR CODE
}

//Hash per tutti 1 su spaced qmer
inline static void GetHash(const string& s_Str, size_t startQmer, size_t length, Hash_Err& hash_err, hash_type (*fConvertion)(char))
{
	hash_err.reset();
//	#pragma omp parallel for ordered
	for(size_t i = startQmer; i < startQmer + length; ++i)
	{
		hash_type ch = (*fConvertion)(s_Str[i]);
//		#pragma omp ordered
		if(ch == 4) //Errore conversione
			hash_err.push_back_error(i);
		else
			hash_err.hash |= ch << ((i - startQmer) * 2);//OR possibile perchè sommo potenze di 4, OR su posizioni diverse, non c'è riporto
	}
}

//Hash per spaced qmer con *
inline static void GetHash(const string& s_Str, size_t startQmer, const SpacedQmer& spaced_qmer,
		Hash_Err& hash_err, hash_type (*fConvertion)(char))
{
	hash_err.reset();
	const Position& pos_one = spaced_qmer.GetPosOne();
	for(size_t j = 0; j < pos_one.size(); ++j)
	{
		hash_type ch = (*fConvertion)(s_Str[startQmer+pos_one[j]]);
		if(ch == 4) //Errore conversione
			hash_err.push_back_error(j);
		else
			hash_err.hash |= ch << (j * 2);//OR possibile perchè sommo potenze di 4, OR su posizioni diverse, non c'è riporto
	}
}

//Hash per spaced qmer con *
inline static void GetHashFromPosOne(const string& s_Str, size_t startQmer, const Position& pos_one,
		Hash_Err& hash_err, hash_type (*fConvertion)(char))
{
	hash_err.reset();
	for(size_t j = 0; j < pos_one.size(); ++j)
	{
		hash_type ch = (*fConvertion)(s_Str[startQmer+pos_one[j]]);
		if(ch == 4) //Errore conversione
			hash_err.push_back_error(j);
		else
			hash_err.hash |= ch << (j * 2);//OR possibile perchè sommo potenze di 4, OR su posizioni diverse, non c'è riporto
	}
}


//Hash veloce con spaced qmer tutti 1
inline static void GetHashes_speedup_previous(const string& s_Str, size_t length,
		Hash_Err_V& vHash, hash_type (*fConvertion)(char)) {
	vHash.clear();
	if(s_Str.size() >= length)
	{
		size_t n_hashes = s_Str.size() - length + 1;
		vHash.resize(n_hashes); //Crea vettore

		GetHash(s_Str, 0, length, vHash[0], fConvertion);//primo da computare a parte
		for(size_t pos=1; pos < vHash.size(); ++pos)
		{
			Hash_Err& prev_hash = vHash[pos-1];
			Hash_Err& curr_hash = vHash[pos];

			//copia hash e sottrai una posizione dal precedente
			curr_hash.hash = prev_hash.hash;
			curr_hash.hash >>= 2; //sposta 2 bit, esce una lettera
			curr_hash.sub_pos_err(1, prev_hash);

			hash_type enter = (*fConvertion)(s_Str[pos+length-1]);
			if(enter == 4)
				curr_hash.push_back_error(length-1);
			else
				curr_hash.hash |= enter << ((length - 1) * 2);	//aggiungi ultimo elemento OR possibile perchè prima ho
																//diviso per 4 e la posizione dove scrivo ha sicuramente 0
		}
	}
}

inline static void GetHashes_naive(const string& s_Str, const SpacedQmer& spaced_qmer,
		Hash_Err_V& vHash, hash_type (*fConvertion)(char))
{
//	bool isAllOne = spaced_qmer.GetWeight() == spaced_qmer.GetQ();
//	if(isAllOne)
//		GetHashes_speedup_previous(s_Str, spaced_qmer.GetQ(), vHash, fConvertion);
//	else
//	{
		vHash.clear();
		if(s_Str.size() >= spaced_qmer.GetQ())
		{
			size_t n_hashes = s_Str.size() - spaced_qmer.GetQ() + 1;
			vHash.resize(n_hashes); //Crea vettore
			//#pragma omp parallel for //num_threads(8)
			for(size_t pos=0; pos < vHash.size(); ++pos)
				GetHash(s_Str, pos, spaced_qmer, vHash[pos], fConvertion);
		}

//	}
}

inline static void compute_hash_for_speedup_previous(const string& s_Str,
													const Position& pos_one_current, const Position& pos_one_prev,
													const PreviousShift& curr_sp_shift,
													const Hash_Err& prev_hash_err,
													size_t idx_curr_hash, Hash_Err& curr_hash_err,
													hash_type (*fConvertion)(char))
{
	curr_hash_err.hash = prev_hash_err.hash; //Copia hash
	curr_hash_err.hash >>= 2*curr_sp_shift.one_exit;//Shifta correttamente
	if(!curr_sp_shift.one_to_remove.empty())
	{
		hash_type reset_one = 0;
		for(size_t j = 0; j < curr_sp_shift.one_to_remove.size(); ++j)
			reset_one |= (hash_type)3 << (curr_sp_shift.one_to_remove[j] * 2);
		curr_hash_err.hash &= ~reset_one;
	}
	if(!prev_hash_err.isCorrect())
	{
		long curr_pos_one = 0;
		for(size_t e = 0; e < prev_hash_err.size_error(); ++e)
			if((curr_pos_one = prev_hash_err[e]-curr_sp_shift.one_exit) >= 0)
				if(pos_one_prev[prev_hash_err[e]]-curr_sp_shift.shift_min == pos_one_current[curr_pos_one])
					curr_hash_err.push_back_error(curr_pos_one);//aggiorna posizione errore
	}


	//aggiorna posizioni da cambiare su hash
	for(size_t j = 0; j < curr_sp_shift.one_to_change.size(); ++j)
	{
		const size_t& i_to_change = curr_sp_shift.one_to_change[j];
		size_t index_char = idx_curr_hash+pos_one_current[i_to_change];
		// per chiascuno di quelli da cambiare faccio la nuova conversione
		hash_type ch = (*fConvertion)(s_Str[index_char]);
		if(ch == 4) //Errore conversione
			curr_hash_err.push_back_error(i_to_change);
		else
			curr_hash_err.hash |= ch << (i_to_change * 2);//OR possibile perchè sommo potenze di 4, OR su posizioni diverse, non c'è riporto
	}

	//aggiorna rimanenti posizioni da cambiare su hash (quelle uscite son già rimosse)
	//TODO:	si elimina questo pezzo (salta if) se il numero di uno son diversi,
	//in quanto non so dove devo andar ad inserire e rimuovere,
	//NB: l'informazione dove inserire e rimuovere è contenuta tutta
	//sui vettori one_to_change e one_to_remove in quest'ultimo caso


	// aggiunge valori alla estremità destra dell seed(in realtà a sinistra del hash)
	// quelli che non potrei calcolare con sovrapposizioni
	if(pos_one_current.size() == pos_one_prev.size())
			for(size_t j = pos_one_current.size()-curr_sp_shift.one_exit; j < pos_one_current.size(); ++j)
			{
				size_t index_char = idx_curr_hash+pos_one_current[j];
				hash_type ch = (*fConvertion)(s_Str[index_char]);
				if(ch == 4) //Errore conversione
					curr_hash_err.push_back_error(j);
				else
					curr_hash_err.hash |= ch << (j * 2);//OR possibile perchè sommo potenze di 4, OR su posizioni diverse, non c'è riporto
			}
	////////////////////////////////////////////////////////////////////
	if(!curr_hash_err.isCorrect())
		curr_hash_err.sort_uniq_err();
}


inline static void GetHashes_speedup_previous(const string& s_Str, const SpacedQmer& spaced_qmer,
		Hash_Err_V& vHash, hash_type (*fConvertion)(char)) {
//	bool isAllOne = spaced_qmer.GetWeight() == spaced_qmer.GetQ();
//	if(isAllOne)
//		GetHashes_speedup_previous(s_Str, spaced_qmer.GetQ(), vHash, fConvertion);
//	else
//	{
		auto get_hash = [&](size_t curr_idx_hash, const PreviousShift& curr_shift){
			Hash_Err& curr_hash = vHash[curr_idx_hash];

			if(spaced_qmer.GetWeight() < curr_shift.GetSize())
				GetHash(s_Str, curr_idx_hash, spaced_qmer, curr_hash, fConvertion);
			else
			{
				size_t pos_hash_get = curr_idx_hash-curr_shift.shift_min;
				const Hash_Err& prev_hash = vHash[pos_hash_get];
				compute_hash_for_speedup_previous(s_Str,
						spaced_qmer.GetPosOne(), spaced_qmer.GetPosOne(),
						curr_shift,
						prev_hash,
						curr_idx_hash, curr_hash,
						fConvertion);
			}
		};
		long n_hashes = s_Str.size()-spaced_qmer.GetQ()+1;
		vHash.clear();
		if(n_hashes>0)
		{
			const V_PreviusShift& shift = spaced_qmer.GetShiftMinChange();
			vHash.resize(n_hashes); //Crea vettore

			GetHash(s_Str, 0, spaced_qmer, vHash[0], fConvertion); // primo da computare a parte
			size_t lim_max = vHash.size();
			size_t lim_min = shift.size() < lim_max ? shift.size() : lim_max;
			for(size_t i = 1; i < lim_min; ++i)//Per tutte le posizioni che contemplano gli shift nel primo pezzo di sequenza
				{
					get_hash(i, shift[i]);
				}
			for(size_t i = lim_min; i < lim_max; ++i)
			{
				get_hash(i, shift.back());
			}
		}
//	}
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Funzione ispirata alla lambda get_hash in GetHashes_speedup_previous che però
// a sua volta non chiama nessuna altra funzione. Come prima cosa scorre il
// gruppo di shift precedenti e da ciascuno recupera le posizioni
// precedentemente identificate.
// Successivamente riempie i buchi rimanenti calcolando la corrispondenza delle
// posizioni rimaste tramite la funzione di codifica.

// A questa funzione vengono passati:

// curr_idx_hash indice dell'hash che voglio calcolare

// shifts gruppo di informazioni sugli hash precedenti che mi permettono di
// recuperare le posizioni di interesse da ciascun hash.

// vHash vettore degli hash dove alla posizione curr_idx_hash andrò a scrivere
// il valore dell'hash corrente che sto calcolando e da dove recupererò i valori
// degli hash precedenti.

// pos_one vettore delle posizioni del seed in cui ho gli uno, è necessario per
// gestire gli errori presenti nella sequenza letta(errori non testati).
// Va utilizzato in abbinata a one_to_keep di ciascun hash precedente che si
// considera

// s_Str sequenza di DNA di cui devo fare l'hash

// fConvertion puntatore alla funzione di conversione da basi azotate a numeri

inline static void compute_hash_with_ISSH(											size_t curr_idx_hash,
																					const V_PreviusShift& shifts,
																					Hash_Err_V& vHash,
																					const Position& pos_one,
																					const string& s_Str,
																					hash_type (*fConvertion)(char))
{
	// creo un puntatore all'hash che voglio calcolare al momento
	// in modo da non passare ogni volta per il vettore
	Hash_Err& curr_hash = vHash[curr_idx_hash];

	// Inizializzo il valore dell'hash che devo calcolare a zero
	curr_hash.hash = 0;
	// Vado a prendere le posizioni che non posso recuperare tramite gli hash
	// precedenti e che andrò a calcolare tramite la funzione di hash alla fine.
	// Questo vettore sicuramente conterrà l'ultima posizione.
	// È stato inserito dentro a shift[0].one_to_change perché era un vettore di posizioni
	// presente in FSH e che non veniva utilizzato in ISSH
	const Position& pos_not_covered_yet = shifts[0].one_to_change;
	// faccio un ciclo per ciascuno degli hash del gruppo da cui devo
	// recuperare posizioni
	for(size_t k = 0; k < shifts.size(); k++)
		{
		//identifico indice dell'hash precedente che mi serve
		size_t pos_hash_get = curr_idx_hash-shifts[k].shift_min;
		//carico l'hash calcolato precedentemente
		const Hash_Err& prev_hash = vHash[pos_hash_get];
		// in partial_hash salvo l'hash precedente da cui voglio recuperare posizioni
		hash_type partial_hash = prev_hash.hash;
		const PreviousShift& curr_sp_shift= shifts[k];
		// faccio lo shift dell'hash in modo da avere la giusta sovrapposizione
		partial_hash >>= 2*curr_sp_shift.one_exit;

		// facendo l'and con la maschera cancello i bit dell'hash che sono
		// sbagliati e mantengo quelli giusti
		partial_hash &= curr_sp_shift.mask;

		// Controllo se l'hash da cui ho preso i valori conteneva errori:
		// con isCorrect torna la lunghezza del vettore di posizioni di errore.
		// Faccio un ciclo fino alla lunghezza tornata da size_error(), prev_hash è
		// di tipo hash_err, al suo interno ha un err_pos e così vado a prenderne la
		// lunghezza. Controllo poi se l'errore del precedente si è propagato sull'
		// hash attuale, se sì inserisco l'errore nella corretta posizione anche in
		// questo hash
		// in realtà se l'input è corretto qui dentro non ci vado mai
		// c'è da risolvere un piccolo problema, va considerato anche il fatto del
		// one_to_keep forse. È comunque da ricontrollare...
		// Questo pseudo controllo è stato inserito perché presente anche nelle
		// altre funzioni di hash. In questo modo viene effettuato sempre lo
		// stesso controllo in ciascuna funzione per non creare differenze di
		// questo tipo, ma se la sequenza in ingresso non contiene errori si
		// potrebbe anche togliere.
		if(!prev_hash.isCorrect())
		{
			long curr_pos_one = 0;
			for(size_t e = 0; e < prev_hash.size_error(); ++e)
				if((curr_pos_one = prev_hash[e]-curr_sp_shift.one_exit) >= 0)
					if(pos_one[prev_hash[e]]-curr_sp_shift.shift_min == pos_one[curr_pos_one])
						curr_hash.push_back_error(curr_pos_one);//aggiorna posizione errore
		}
		curr_hash.hash |= partial_hash;
	}

// A questo punto si devono inserire, calcolati uno ad uno con la funzione di
// codifica, tutti i valori contenuti nelle posizioni di pos_not_covered_yet.
// Dentro questo vettore ho anche la posizione dell'ultimo valore dell'hash
// che va calcolato in ogni caso.
	for(size_t j = 0; j < pos_not_covered_yet.size(); ++j)
	{
		const size_t& i_to_change = pos_not_covered_yet[j];
		// identifico a che indice della sequenza devo calcolare la base
		size_t index_char = curr_idx_hash+pos_one[i_to_change];
		// per ciascuno di quelli da cambiare faccio la nuova conversione tramite la
		// funzione di hash
		hash_type ch = (*fConvertion)(s_Str[index_char]);
		if(ch == 4) //Errore conversione
			curr_hash.push_back_error(i_to_change);
		else
			curr_hash.hash |= ch << (i_to_change * 2);//OR possibile perchè sommo potenze di 4, OR su posizioni diverse, non c'è riporto
	}
};


inline static void GetHashes_with_ISSH(const string& s_Str,
																			 const SpacedQmer& spaced_qmer,
																			 Hash_Err_V& vHash,
																			 hash_type (*fConvertion)(char))
{
  // calcolo il numero totale di hash che calcolerò per questa stinga
	long n_hashes = s_Str.size()-spaced_qmer.GetQ()+1;
	vHash.clear();
	// Comincio a fare i calcoli solo se riesco a fare almeno un hash
	if(n_hashes>0)
	{
		// carico il vettore dei gruppi di shift precedenti che mi serve per
		// gestire sia il caso a regime sia tutto il transitorio.
		const V_V_PreviusShift& V_shifts = spaced_qmer.GetMultipleShifts();
		vHash.resize(n_hashes);
		// per lo 0-esimo hash ovviamente non posso recuperare niente, lo calcolo
		// posizione per posizione.
		GetHash(s_Str, 0, spaced_qmer, vHash[0], fConvertion);
		// Il transitorio mi basta farlo dal primo alla lunghezza del vettore di vettori
		// di shift che voglio andare a utilizzare che gestisce il transitorio dall'1
		// fino alla fine. Ad ogni passo del trasitorio chiamo la stessa funzione
		// passandogli il gruppo di hash precedenti corrispondente
		const Position& pos_one = spaced_qmer.GetPosOne();
		size_t i=1;
		//transitorio
		for(; i<(size_t)n_hashes && i<V_shifts.size(); i++)
		{
			// se il vettore di shift è vuoto è sicuro che non posso recuperare nessuna
			// posizione. Per evitare controlli e passaggi inutili si utilizza la funzione
			// posizione per posizione standard (gestione migliore per il caso 10101...)
			if(V_shifts[i].size()==0)
			{
				GetHash(s_Str, i, spaced_qmer, vHash[i], fConvertion);
			}
			else
			{
				compute_hash_with_ISSH(i, V_shifts[i], vHash, pos_one, s_Str, fConvertion);
			}
		}

		// regime nel caso a regime si passa alla funzione compute sempre lo stesso
		// gruppo di hash
		for(; i < (size_t)n_hashes; i++)
		{
			// in questo caso si suppone che a regime siano recuperabili tutte le posizioni
			// quindi non ha senso vedere se il vettore di posizioni precedenti è
			// vuoto
			compute_hash_with_ISSH(i, V_shifts[0], vHash, pos_one, s_Str, fConvertion);
		}
	}
}


#endif /* HASHFUNCTION_H_ */
