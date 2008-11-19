#ifndef IDXVECTOR_H_
#define IDXVECTOR_H_

#include <algorithm>
#include <deque>

namespace SLOM {



// An indexed vector
//TODO override erase() etc.
template<typename T>
class IdxVector : public std::deque<T*>{
	typedef std::deque<T*> Container;
	/** helper function for find */
	static inline bool compareIdx(T* x, int idx){
		return x->idx < idx;
	}
	int lastIdx;
	
public:
	IdxVector() : lastIdx(0) {}
	void push_back(T *m){
		m->idx=lastIdx;
		Container::push_back(m);
		lastIdx += m->getDim();
	}
	int getDim() const { return lastIdx; }
	
	
	
	/** 
	 * find finds the entry, which starts at index idx.
	 * TODO this isn't used.
	 */ 
	std::pair<T*, int> find(int idx) const {
		T* x = *std::lower_bound(this->begin(), this->end(), idx, compareIdx );
		return std::make_pair(x, idx - x->idx);
	}
};


}  // namespace SLSQ

#endif /*IDXVECTOR_H_*/
