#ifndef __CORRELATOR_MULTITAU_H_INCLUDED__
#define __CORRELATOR_MULTITAU_H_INCLUDED__

/**
 * Correlator class implementing a multi-tau correlator.
 * This class is NOT intended to be optimized or fast,
 * it tries to build a hardware equivalent of a parallel
 * multi-tau correlator without beeing parallel. ;-)
 * This is how it should look like (for a (3,2) invocation):
 *
 *
 */
class block;
class lag;

class correlatorjb
{
	private:
		block **blocks;/**Blocks running at different frquencies*/
		unsigned int block_count;/**Number of blocks in current correlator setup*/
		correlatorjb *next;/**One way links list, reference to the next block in chain.*/
	public:
		correlatorjb(unsigned int,unsigned int);
		~correlatorjb();
		void run(int,int);
		void print(unsigned int*);
		unsigned int		get_channel_count();
		int *get_array_lin();
		float *get_array_lin_G();
		int **get_array();
		float **get_array_G();
		unsigned long long sum;
};

/**
 * Class block, implements one block of a multi-tau correlator.
 * One blocks runs at the same clock.
 */
class block
{
	private:
		lag **lags;
		unsigned int lag_count;
		block *next;
		/**Number of time steps processed*/
		unsigned int ticks;
		/**Id of the current block. Used for printing purpose.*/
		unsigned int id;
	public:
		block(unsigned int,unsigned int);
		~block();
		void run(int,int);
		void print(unsigned int *);
		void set_next(block*);
		int get_lag_count();
		void fill_array_lin(int *array);
		int		get_sum(unsigned int lag_number);
		void fill_array(int **,unsigned int, unsigned int *);
		unsigned int get_ticks() { return ticks; }
};

struct values
{
	int cur,old;
};

class lag {
	private:
		/**This is used to implemnt the flip flop.*/
		values				value;
		/**Next lag in chain.*/
		lag						*next;
		/**Id of this lag in curent block.*/
		unsigned int 	id;
  public:
		/**Result of correlation*/
		int sum;
		/**Sum of the last two time steps. Used for connection to next block.*/
		int acc;

				lag(unsigned int,void *);
	void	run (int,int);
	void	set_next (lag *);
	void	print(unsigned int);
	int		get_sum();
};


#endif
