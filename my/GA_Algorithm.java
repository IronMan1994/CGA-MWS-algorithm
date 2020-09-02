package my;

import java.lang.reflect.Method;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

public class GA_Algorithm
{
	//Data matrix
	public List<Data_Array> data_Array;
	//Number of genes in the data
	int n;
	//Gene set size
	int k;	
	//Population size
	int geneSize;
	//Total iterations
	int genestep;
	// Number of times to run the algorithm
	int MaxStep;
	//Control iteration steps
	int Step;
	//Mutation probability
	double pm;
	//Gene name
	public String[] name;
	//Initial sequence of genes
	public int[] name_index;
	//The optimal value obtained by running the algorithm several times
	public List<SpeciesIndividual_index> maxgene;
	//Contemporary population
	public List<SpeciesIndividual> pop;
	//Next generation population
	public List<SpeciesIndividual> nextPOP;
	//Global random variables
	public Random random;
	//Global computing model
	public Method method;
	//time
	public long start;
	public long end;
	//Tools
	public Tools tools;
	
	
	//Initialization data, variables
	public void initData(String[] paths, int n, int k, int geneSize, int genestep, double pm, int step, String methodName, Tools tools) throws Exception
	{
		this.tools = tools;
		//初始化各种变量
		this.n = n;
		this.k = k;
		this.geneSize = geneSize;
		this.genestep = genestep;
		this.MaxStep = step;
		this.Step = 0;
		this.pm = pm;
		this.name = new String[this.n];
		this.name_index = new int[this.n];
		this.maxgene = new ArrayList<SpeciesIndividual_index>();
		this.pop = new ArrayList<SpeciesIndividual>();
		this.nextPOP = new ArrayList<SpeciesIndividual>();
		this.random = new Random(System.currentTimeMillis());	
		this.data_Array = new ArrayList<Data_Array>();
		
		tools.initData(paths, data_Array, n, name, name_index);
		
		//Global computing model
		Method[] methods = this.getClass().getMethods();
		for(Method m : methods)
		{
			if(methodName.equals(m.getName()))
			{
				this.method = m;
				break;
			}
		}
	}
	
	
	//Constructing initial population
	public void createBeginningSpecies() throws Exception
	{
		pop.clear();
		int count = 0;
		List<Integer> index = new ArrayList<Integer>();
		for(int i = 0; i < name_index.length; i++){
			index.add(i);
		}
		while(count < this.geneSize)
		{
			Collections.shuffle(index);
			int[] chromosome = new int[k];
			for(int m = 0; m < k; m++)
			{
				chromosome[m] = index.get(m);
			}
			
			double[] fitness_result;
			Object[] args = {chromosome, this.data_Array.get(0)};
			fitness_result = (double[]) method.invoke(this, args);
			SpeciesIndividual speciesIndividual = new SpeciesIndividual(chromosome, fitness_result);
			pop.add(speciesIndividual);
			count++;
		}//while	
	}
	
	//Selecting individuals for iterative operations
	public void select()
	{
		nextPOP.clear();
		tools.sort_pop_list(pop);
		//最好的
		for(int i = 0; i < geneSize; i++){
			nextPOP.add(pop.get(i));
		}
		tools.copy_list(nextPOP, pop);
	}
	
	
	//Calculate roulette selection probability
	public double[] crossover_rate(List<SpeciesIndividual> Cross_pool)
	{
		int geneSize = Cross_pool.size();
		int fuzhiNum = 0;
		//Roulette probability array
		double[] rate_array = new double[geneSize - fuzhiNum + 1];
		double rate = 0.0;
		double sum = tools.sum_fitness(fuzhiNum, geneSize, Cross_pool, geneSize);
		rate_array[0] = 0.0;
		for(int i = 1; i < rate_array.length; i++)
		{
			rate += Cross_pool.get(fuzhiNum - 1 + i).fitness[0] / sum;
			rate_array[i] = rate;
		}
		return rate_array;
	}
	
	//Crossover operation
	public void crossover(List<SpeciesIndividual> Cross_pool, double[] rate_array) throws Exception
	{
		//Selecting Parental Individuals
		double r1 = Math.random();
		int index1 = 0;
		SpeciesIndividual parent1 = null;
		for(int j = 0; j < rate_array.length - 1; j++)
		{
			if(r1 >= rate_array[j] && r1 <= rate_array[j + 1])
			{
				parent1 = Cross_pool.get(j);
				index1 = j;
			}
		}
		
		SpeciesIndividual parent2 = null;
		int index2 = index1;
		while(index1 == index2)
		{
			double r2 = Math.random();
			for(int j = 0; j < rate_array.length - 1; j++)
			{
				if(r2 >= rate_array[j] && r2 <= rate_array[j + 1])
				{
					parent2 = Cross_pool.get(j);
					index2 = j;
				}
			}
		}
		
		//Start the crossover operation
		int[] chromosome1 = new int[this.k];
		int[] chromosome2 = new int[this.k];
		List<Integer> remainList1 = new ArrayList<Integer>();
		List<Integer> remainList2 = new ArrayList<Integer>();	
		
		Set<Integer> temp_set = new HashSet<Integer>();
		for(int i = 0; i < parent1.chromosome.length; i++){
			temp_set.add(parent1.chromosome[i]);
			remainList1.add(parent1.chromosome[i]);
		}
		
		int count = 0;
		for(int i = 0; i < parent2.chromosome.length; i++){
			try{
				if(!temp_set.add(parent2.chromosome[i])){
					remainList1.remove(remainList1.indexOf(parent2.chromosome[i]));
					chromosome1[count] = parent2.chromosome[i];
					chromosome2[count] = parent2.chromosome[i];
					count++;
				}
				else{
					remainList2.add(parent2.chromosome[i]);
				}
			}
			catch(Exception e)
			{
				e.printStackTrace();
			}
		}
			
		List<Integer> temp_list = new ArrayList<Integer>();	
		for(Integer a : remainList1){
			temp_list.add(a);
		}
		for(Integer a : remainList2){
			temp_list.add(a);
		}
		Collections.shuffle(temp_list);
		for(int i = 0; i < temp_list.size(); i+=2)
		{
			boolean wi = random.nextBoolean();
			if(wi == false){
				chromosome1[count] = temp_list.get(i);
				chromosome2[count] = temp_list.get(i+1);
			}
			else{
				chromosome2[count] = temp_list.get(i);
				chromosome1[count] = temp_list.get(i+1);
			}
			count++;
		}
		Object[] args = {chromosome1, this.data_Array.get(0)};
		double[] fitness_result = (double[]) method.invoke(this, args);
		SpeciesIndividual speciesIndividual = new SpeciesIndividual(chromosome1, fitness_result);
		
		Object[] args2 = {chromosome2, this.data_Array.get(0)};
		double[] fitness_result2 = (double[]) method.invoke(this, args2);
		SpeciesIndividual speciesIndividual2 = new SpeciesIndividual(chromosome2, fitness_result2);
	
		pop.add(speciesIndividual);
		pop.add(speciesIndividual2);
	}
	
	//Mutation operation
	public void mutate_SA(int index) throws Exception
	{
		double rate = Math.random();
		if(rate <= this.pm){
			mutate(index);
		}
	}
	
	//Mutation operation
	public void mutate(int Spindex) throws Exception
	{
		int[] chromosome = new int[k];
		chromosome = pop.get(Spindex).chromosome;
		double init_fitness = pop.get(Spindex).fitness[0];
		
		//Randomly delete a gene
		int[] max_gene = new int[this.k];
		max_gene[this.k-1] = -1;
		int removeIndex = random.nextInt(65535) % k;
		int count2 = 0;
		for(int l = 0; l < chromosome.length; l++)
		{
			if(l != removeIndex){
				max_gene[count2++] = chromosome[l];
			}
		}
		//Searching for candidate gene set
		List<Integer> houxuan = new ArrayList<Integer>();
		int count = 0;
		boolean ifcunzai;
		for(int j = 0; j < name_index.length; j++){	
			ifcunzai = false;
			for(int k = 0; k < chromosome.length; k++){
				if(name_index[j] == chromosome[k]){
					ifcunzai = true;
					break;
				}
			}
			if(ifcunzai == false){
				houxuan.add(name_index[j]);
				count++;
			}			
		}
		double maxfitness = -Double.MAX_VALUE;
		double tempfitness = 0.0;
		int maxIndex = 0;
		double[] max_result = null;
		Collections.shuffle(houxuan);
		int random_count = (int) Math.sqrt(houxuan.size());
		int step = 0;
		while(step < random_count){
			max_gene[this.k-1] = houxuan.get(step);
			Object[] args1 = {max_gene, this.data_Array.get(0)};
			double[] tempfitness_result = (double[]) method.invoke(this, args1);
			tempfitness = tempfitness_result[0];
			if(maxfitness < tempfitness){
				maxfitness = tempfitness;
				maxIndex = step;
				max_result = tempfitness_result;
			}
			step++;
		}
		//Decide whether to replace the individual
		if(maxfitness > init_fitness){
			max_gene[this.k - 1] = houxuan.get(maxIndex);
			pop.get(Spindex).chromosome = max_gene;
			pop.get(Spindex).fitness = max_result;
		}
	}
		
	
	//Fitness calculation based on coefficient of variation
	public double[] calfitness_Cov(int[] chromosome, Data_Array B)
	{	
		double[] result = new double[5];
		double fitness = 0.0;
		double red = 0.0,ifcancer = 0.0;
		int x = 0;
		int y = 0;
		List<float[]> sample_one = B.A;	
		for(int i = 0; i < sample_one.size();i++)
		{
			float[] temp = B.A.get(i);
			double[] temp_mean = new double[chromosome.length];
			double sum_line = 0.0; 
			double max = 0;
			for(int j = 0; j < chromosome.length;j++)
			{
				if(temp[chromosome[j]] > max)
				{
					max = temp[chromosome[j]];
				}	
				if(temp[chromosome[j]] >= 1){
					y++;
				}
				sum_line += temp[chromosome[j]];
				temp_mean[j] = temp[chromosome[j]];
			}
			if(max >= 1){
				x++;
			}
			if(sum_line != 0){
				double mean = sum_line / chromosome.length;
				double fangcha_fenzi = 0.0;
				for(int j = 0; j < temp_mean.length; j++){
					fangcha_fenzi += Math.pow(temp_mean[j] - mean, 2);
				}
				double cov = Math.sqrt(fangcha_fenzi / (this.k-1)) / mean;
				
				if(max <= 0.5){
					cov = cov / (Math.sqrt(this.k)*2) ;
				}
				red += cov;
			}
			ifcancer += max;	
		}
		BigDecimal bd = new BigDecimal(ifcancer).setScale(4,RoundingMode.UP);
		ifcancer = bd.doubleValue();
		
		BigDecimal bd2 = new BigDecimal(red).setScale(4,RoundingMode.UP);
		red = bd2.doubleValue();
		fitness = (ifcancer + red);
		
		BigDecimal bd3 = new BigDecimal(fitness).setScale(4,RoundingMode.UP);
		fitness = bd3.doubleValue();
		
		BigDecimal bd4 = new BigDecimal(ifcancer).setScale(4,RoundingMode.UP);
		
		result[0] = fitness;
		result[1] = x;
		result[2] = y;
		result[3] = bd4.doubleValue();
		result[4] = Double.valueOf(red);
		
		return result;	
	}
	

	//Based on the original maximum weight submatrix model fitness calculation.
	public double[] calfitness_01(int[] chromosome, Data_Array B)
	{
		double[] result = new double[3];
		double fitness = 0.0;
		double red = 0.0,ifcancer = 0.0, overlap = 0.0;
		List<float[]> sample_one = B.A;	
		for(int i = 0; i < sample_one.size();i++)
		{
			float[] temp = B.A.get(i);
			int tempifcancer = 0;
			for(int j = 0; j < chromosome.length;j++)
			{
				if(temp[chromosome[j]] != 0)
				{
					red += temp[chromosome[j]];
					tempifcancer = 1;
				}
			}
			ifcancer += tempifcancer;		
		}
		fitness = 2 * ifcancer - red;
		overlap = red - ifcancer;

		result[0] = fitness;
		result[1] = Double.valueOf(ifcancer);
		result[2] = Double.valueOf(red);
		
		return result;
	}	
	
}










