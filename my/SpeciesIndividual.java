package my;


public class SpeciesIndividual
{
	//一条染色体
	public int[] chromosome;

	

	//几个对应的适应度(第一个为总和的适应度)
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
