/*
 * Declarations for Deme class to evolve a genetic algorithm for the
 * travelling-salesperson problem.  A deme is a population of individuals.
 */

#include "chromosome.hh"
#include "deme.hh"
#include <random>

// Generate a Deme of the specified size with all-random chromosomes.
// Also receives a mutation rate in the range [0-1].
Deme::Deme(const Cities* cities_ptr, unsigned pop_size, double mut_rate)
{
  std::vector<Chromosome*> pop_(pop_size);


  for (int i = 0; i < pop_size; i++){
  	//creates a chromosome based on cities_ptr.
  	Chromosome* newGuy = new Chromosome(cities_ptr);

  	pop_[i] = newGuy;
  }
  mut_rate_ = mut_rate;



}

// Clean up as necessary
Deme::~Deme()
{
  for(int i = 0; i < pop_.size(); i++){
  	delete pop_[i];
  }
  delete this;
}

// Evolve a single generation of new chromosomes, as follows:
// We select pop_size/2 pairs of chromosomes (using the select() method below).
// Each chromosome in the pair can be randomly selected for mutation, with
// probability mut_rate, in which case it calls the chromosome mutate() method.
// Then, the pair is recombined once (using the recombine() method) to generate
// a new pair of chromosomes, which are stored in the Deme.
// After we've generated pop_size new chromosomes, we delete all the old ones.
/*void Deme::compute_next_generation()
{
  // find 2 parents
  Chromosome* parent1 = select_parent();
  Chromosome* parent2 = select_parent();

  // take random number and check for mutation
  int mut_num = rand() % (100*1000); // 100,000
  									 // can make that line a lambda function...
  double mut1 = mut_num / 100000.0;

  int mut_num = rand() % (100*1000);
  double mut2 = mut_num / 100000.0;  
  
  if(mut1<mut_rate) {
  	parent1->mutate();
  }
  if(mut2<mut_rate) {
  	parent2->mutate();
  }

  std::pair<Chromosome*, Chromosome*> children;
  children = parent1 -> recombine(parent2);

  //put children into vector and delete old parents

}*/

void Deme::compute_next_generation()
{
  // make 2 parents
  for (int i = 0; i < pop_.size(); i++){ 
    Chromosome* parent1 = select_parent();
    Chromosome* parent2 = select_parent();
    int rand_val1 = rand() % (100*1000); 
    double mut1 = mut_num / 100000.0;
    int mut_num = rand() % (100*1000);
    double rand_val2 = mut_num / 100000.0;  

    if (rand_val1 < mut_rate_){
      parent1->Chromosome::mutate();
    }
    if (rand_val2 < mut_rate_){
      parent1->Chromosome::mutate();
    }
    std::pair<Chromosome*, Chromosome*> new_children = parent1->recombine(parent2);
  }

  // take random number and check for mutation
  
}


// Return a copy of the chromosome with the highest fitness.
const Chromosome* Deme::get_best() const
{
  Chromosome* mostFit = pop_[0];
  double shortestDist = mostFit->get_fitness();

  for(int i = 1; i < pop_.size(); i++){
  
  	if (pop_[i]->get_fitness()>shortestDist){
  		mostFit = pop_[i];
  	}
  
  }
  return mostFit;
}

// Randomly select a chromosome in the population based on fitness and
// return a pointer to that chromosome.
Chromosome* Deme::select_parent()
{

  double halfwayPt = get_best()->get_fitness()/2;

  // if first random guy isn't in the upper 50%, we try again for another
  // fit guy, and just select the next random guy as our parent
  int rando = rand() % pop_.size();
  return pop_[rando > halfwayPt ? rando : rand() % pop_.size()];
}
