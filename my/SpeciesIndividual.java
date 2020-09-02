package my;


public class SpeciesIndividual
{
	//chromosome
	public int[] chromosome;

	//fitness(the first is the total fitness)
	public double[] fitness;
	
	public SpeciesIndividual(int[] chromosome, double[] fitness)
	{
		this.chromosome = chromosome;
		this.fitness = fitness;
	}
	
	public SpeciesIndividual()
	{
		
	}

	@Override
	public String toString()
	{
		String str = "";
		for(int i = 0; i < chromosome.length; i++)
		{
			str += chromosome[i] + "\t";
		}
		for(int i = 0; i < fitness.length; i++)
		{
			str += fitness[i] + "\t";
		}
		return str;
	}	
}
