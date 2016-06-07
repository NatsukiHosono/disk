#pragma once

template <class TreePtcl> class Tree{
	static const std::size_t MAX_LEVEL = 20;
	static const std::size_t N_LEAF = 16;
	static const double EXPAND = 1.0;
	static const double MAC = 0.5;//Multipole Acceptance Criterion
	typedef unsigned long long int ulli;
	struct node_t{
		vec3<double> mass_center, min, max;
		//mat33<double> Quad;
		double mass, size, radius, dmax;
		std::size_t Nptcl, head_ptcl, more;
		unsigned short int Nbranch;
		unsigned int level;
		node_t(){
			Nptcl = 0;
			mass_center = 0;
			mass = 0;
		}
		~node_t(){
		}
		inline bool isEmpty() const{
			return (Nptcl == 0) ? true : false;
		}
		inline bool isLeaf(const unsigned int N) const{
			return (Nptcl <= N) ? true : false;
		}
		inline const node_t* const getPointer() const{
			return this;
		}
		inline node_t* getPointer(){
			return this;
		}
	};
	struct ptcl_ptr_t{
		ulli key;
		const TreePtcl* ptr;
		bool operator < (const ptcl_ptr_t& right) const{
			return key < right.key;
		}
		ptcl_ptr_t(const ulli _key, const TreePtcl* const _ptr) : key(_key), ptr(_ptr){
		}
	};
	//member variables;
	std::vector<node_t> node;
	std::vector<ptcl_ptr_t> ptcl_ptr;
	//member functions;
	node_t createRoot(const std::vector<TreePtcl >& ptcl, const kernel_t<double>& kernel){
		node_t root;
		root.min = + 1.0e+30;
		root.max = - 1.0e+30;
		root.mass_center = 0;
		root.mass        = 0;
		root.Nptcl = ptcl.size();
		root.head_ptcl = 0;
		root.level  = 0;
		root.radius = 0;
		root.dmax   = 0;
		for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
			root.min.x = math::min(root.min.x, ptcl[i].r.x - kernel.width * ptcl[i].h);
			root.min.y = math::min(root.min.y, ptcl[i].r.y - kernel.width * ptcl[i].h);
			root.min.z = math::min(root.min.z, ptcl[i].r.z - kernel.width * ptcl[i].h);
			root.max.x = math::max(root.max.x, ptcl[i].r.x + kernel.width * ptcl[i].h);
			root.max.y = math::max(root.max.y, ptcl[i].r.y + kernel.width * ptcl[i].h);
			root.max.z = math::max(root.max.z, ptcl[i].r.z + kernel.width * ptcl[i].h);
		}
		root.size = math::max(math::max((root.max - root.min).y, (root.max - root.min).z), (root.max - root.min).x);
		return root;
	}
	//key maker
	inline ulli SpreadBits(ulli x) const{
		x = (x | (x << 32)) & 0x7fff00000000ffff; // 0b0111111111111111000000000000000000000000000000001111111111111111
		x = (x | (x << 16)) & 0x00ff0000ff0000ff; // 0b0000000011111111000000000000000011111111000000000000000011111111
		x = (x | (x <<  8)) & 0x700f00f00f00f00f; // 0b0111000000001111000000001111000000001111000000001111000000001111
		x = (x | (x <<  4)) & 0x30c30c30c30c30c3; // 0b0011000011000011000011000011000011000011000011000011000011000011
		x = (x | (x <<  2)) & 0x1249249249249249; // 0b0001001001001001001001001001001001001001001001001001001001001001
		return x;
	}
	inline ulli zorder3d(const ulli x, const ulli y, const ulli z) const{
		return (SpreadBits(x) | (SpreadBits(y) << 1) | (SpreadBits(z) << 2));
	}
	inline unsigned int GetBranchId(const ulli key, const unsigned int level) const{
		return key >> (3 * MAX_LEVEL - 3 * level) & 0x7; // 0b111
	}
	inline bool isOverwrapping(const TreePtcl& ptcl, const kernel_t<double>& kernel, const node_t& node, const boundary_t<double>& boundary) const{
		return (abs(boundary.Periodic(ptcl.r - node.mass_center)) > math::max(node.radius, node.dmax + EXPAND * kernel.width * ptcl.h)) ? false : true;
	}
	inline bool isClose(const TreePtcl& ptcl, const node_t& node) const{
		return (abs(ptcl.r - node.mass_center) * MAC < node.size) ? true : false;
	}
	//public member functions
	public:
	Tree(){
		node.clear();
		ptcl_ptr.clear();
	}
	//make tree structure
	void Plant(const std::vector<TreePtcl >& ptcl, const kernel_t<double>& kernel){
		//Clear;
		node.clear();
		ptcl_ptr.clear();
		ptcl_ptr.reserve(ptcl.size());
		//Set root domain;
		node.push_back(createRoot(ptcl, kernel));
		//Create particle pointer
		const ulli grid_size = 1 << MAX_LEVEL;
		for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
			struct{
				ulli x, y, z;
			} grid_address;
			grid_address.x = static_cast<ulli>((ptcl[i].r.x - node[0].min.x) / (node[0].max.x - node[0].min.x) * static_cast<double>(grid_size));
			grid_address.y = static_cast<ulli>((ptcl[i].r.y - node[0].min.y) / (node[0].max.y - node[0].min.y) * static_cast<double>(grid_size));
			grid_address.z = static_cast<ulli>((ptcl[i].r.z - node[0].min.z) / (node[0].max.z - node[0].min.z) * static_cast<double>(grid_size));
			const ulli zkey = zorder3d(grid_address.x, grid_address.y, grid_address.z);
			const TreePtcl* const ptr = ptcl[i].getPointer();
			ptcl_ptr.push_back(ptcl_ptr_t(zkey, ptr));
		}
		std::sort(ptcl_ptr.begin(), ptcl_ptr.end());
		//__gnu_parallel::sort(ptcl_ptr.begin(), ptcl_ptr.end());
		
		//Create tree structure;
		for(std::size_t n = 0 ; n < node.size() ; ++ n){
			if(node[n].isLeaf(N_LEAF) == false){
				node[n].more = node.size();//
				//have a baby nodes
				node_t child[8];
				for(short unsigned int c = 0 ; c < 8 ; ++ c){
					child[c].size    = 0.5 * node[n].size;
					child[c].level   = node[n].level + 1;
					child[c].radius  = 0.0;
					child[c].dmax    = 0.0;
					child[c].Nbranch = 0;
					child[c].Nptcl   = 0;
					child[c].more    = 0;//NULL
				}
				for(std::size_t i = node[n].head_ptcl ; i < node[n].head_ptcl + node[n].Nptcl ; ++ i){
					const std::size_t c = GetBranchId(ptcl_ptr[i].key, node[n].level + 1);
					++ child[c].Nptcl;
				}
				child[0].head_ptcl = node[n].head_ptcl;
				for(std::size_t c = 1 ; c < 8 ; ++ c){
					child[c].head_ptcl = child[c-1].head_ptcl + child[c-1].Nptcl;
				}
				node[n].Nbranch = 0;
				for(std::size_t c = 0 ; c < 8 ; ++ c){
					if(__builtin_expect(child[c].isEmpty() == true, 0)) continue;
					++ node[n].Nbranch;
					node.push_back(child[c]);
				}
			}
		}
	}
	//get cell properties
	void setCellProperties(const kernel_t<double>& kernel){
		//Create bounding box;
		#if 0 //Top Down;
		for(vector<node_t>::iterator cell = node.begin() ; cell != node.end() ; ++ cell){
			cell->mass_center = 0;
			cell->mass = 0;
			cell->radius = 0;
			cell->dmax = 0;
			for(std::size_t i = cell->head_ptcl ; i < cell->head_ptcl + cell->Nptcl ; ++ i){
				cell->mass_center += ptcl_ptr[i].ptr->r * ptcl_ptr[i].ptr->Q.mass;
				cell->mass += ptcl_ptr[i].ptr->Q.mass;
			}
			cell->mass_center /= cell->mass;
			for(std::size_t i = cell->head_ptcl ; i < cell->head_ptcl + cell->Nptcl ; ++ i){
				cell->radius = math::max(cell->radius, abs(cell->mass_center - ptcl_ptr[i].ptr->r) + EXPAND * kernel.width * ptcl_ptr[i].ptr->h);
				cell->dmax   = math::max(cell->dmax  , abs(cell->mass_center - ptcl_ptr[i].ptr->r));
			}
		}
		#else //Bottom Up;
		for(typename std::vector<node_t>::reverse_iterator cell = node.rbegin() ; cell != node.rend() ; ++ cell){
			cell->mass_center = 0;
			cell->mass = 0;
			cell->radius = 0;
			cell->dmax = 0;
			if(__builtin_expect(cell->isLeaf(N_LEAF) == true, 0)){
				for(std::size_t i = cell->head_ptcl ; i < cell->head_ptcl + cell->Nptcl ; ++ i){
					cell->mass_center += ptcl_ptr[i].ptr->r * ptcl_ptr[i].ptr->mass;
					cell->mass += ptcl_ptr[i].ptr->mass;
				}
				cell->mass_center /= cell->mass;
				for(std::size_t i = cell->head_ptcl ; i < cell->head_ptcl + cell->Nptcl ; ++ i){
					const double dis = abs(cell->mass_center - ptcl_ptr[i].ptr->r);
					cell->radius = math::max(cell->radius, dis + EXPAND * kernel.width * ptcl_ptr[i].ptr->h);
					cell->dmax   = math::max(cell->dmax  , dis);
				}
			}else{
				for(std::size_t b = cell->more ; b < cell->more + cell->Nbranch ; ++ b){
					cell->mass_center += node[b].mass * node[b].mass_center;
					cell->mass        += node[b].mass;
				}
				cell->mass_center /= cell->mass;
				for(std::size_t b = cell->more ; b < cell->more + cell->Nbranch ; ++ b){
					const double dis = abs(cell->mass_center - node[b].mass_center);
					cell->radius = math::max(cell->radius, dis + node[b].radius);
					cell->dmax   = math::max(cell->dmax  , dis + node[b].dmax  );
				}
			}
		}
		#endif
	}
	//set Ngb List
	void setNeighbourList(std::vector<TreePtcl >& ptcl, const kernel_t<double>& kernel, const boundary_t<double>& boundary){
		#pragma omp parallel for //schedule(guided)
		for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
			ptcl[i].ngb.clear();
			std::stack<const node_t*> stack;
			stack.push(node[0].getPointer());
			while(stack.empty() == false){
				const node_t* const cell = stack.top(); stack.pop();
				if(isOverwrapping(ptcl[i], kernel, *cell, boundary) == true){
					if((*cell).isLeaf(N_LEAF) == true){
						for(std::size_t j = (*cell).head_ptcl ; j < (*cell).head_ptcl + (*cell).Nptcl ; ++ j){
							if(abs(boundary.Periodic(ptcl[i].r - ptcl_ptr[j].ptr->r)) < kernel.width * std::max(ptcl[i].h, ptcl_ptr[j].ptr->h)) ptcl[i].ngb.push_back(ptcl_ptr[j].ptr);
							//ptcl[i].ngb.push_back(ptcl_ptr[j].ptr);
						}
					}else{
						for(std::size_t b = (*cell).more ; b < (*cell).more + (*cell).Nbranch ; ++ b) stack.push(node[b].getPointer());
					}
				}
			}
		}
	}
	//get gravity
	void setSelfGravity(std::vector<TreePtcl >& ptcl);
	//Dumping function
	void Dump(char const * const filename) const{
		FILE* fp = fopen(filename, "w");
		for(std::size_t i = 0 ; i < node.size() ; ++ i){
			fprintf(fp, "%e\t%e\t%e\t%e\t%e\n", node[i].mass_center.x, node[i].mass_center.y, node[i].mass_center.z, node[i].radius, node[i].dmax);
		}
		fclose(fp);
	}
};

