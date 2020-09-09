package my;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class Run
{
	//Cross pool
	public List<SpeciesIndividual> Cross_pool = new ArrayList<SpeciesIndividual>();
	//time
	public long start;
	public long end;
	public Tools tools;
	
	//Update cross pool
	public void update_Cross_pool()
	{
		List<SpeciesIndividual> tempList = new ArrayList<SpeciesIndividual>();
		tools.sort_pop_list(Cross_pool);
		for(int i = 0; i < Cross_pool.size() / 2; i++){
			SpeciesIndividual temp = tools.copy_SpeciesIndividual(Cross_pool.get(i));
			tempList.add(temp);
		}
		tools.copy_add_list(tempList, Cross_pool);
	}
	
	//Execution procedure
	public void run(String[] paths, int n ,int k, int geneSize, int genestep, double pm,int step, int P_step, String methodName)
	{	
		tools = new Tools();
		GA_Algorithm ga1 = new GA_Algorithm();
		GA_Algorithm ga2 = new GA_Algorithm();
		try
		{
			//Initialization data, variables
			//Two threads---Two populations
			ga1.initData(paths, n, k, geneSize/2, genestep, pm, step, methodName, tools);
			ga2.initData(paths, n, k, geneSize/2, genestep, pm, step, methodName, tools);
			List results = new ArrayList<>();
			double one_time = 0.0;
			double time_count = 0.0;
			for(int i = 0; i < ga1.MaxStep; i++)
			{
				start = System.currentTimeMillis();	
				System.out.println(String.format("NO.%-2d Execute the algorithm to reconstruct the initial population", i + 1));
				//threads
				ExecutorService es = Executors.newFixedThreadPool(2);
				Callable<Boolean> callable1 = new Callable<Boolean>() {
					@Override
					public Boolean call() throws Exception
					{
						try{
							int j = 0;
							double[] rate_array = ga1.crossover_rate(Cross_pool);
							double[] aa = ga1.crossover_rate(ga1.pop);
							while(j < ga1.geneSize / 4){
								ga1.crossover(ga1.pop, aa);
								ga1.crossover(Cross_pool, rate_array);
								j++;
							}
							j=0;
							while(j < ga1.geneSize/2){
								ga1.mutate_SA(ga1.geneSize + j);
								ga1.mutate_SA(ga1.geneSize + j + 1);
								j += 2;
							}
							ga1.select();
							return true;
						}
						catch(Exception e){
							e.printStackTrace();
							return false;
						}	
					}
				};		
				Callable<Boolean> callable2 = new Callable<Boolean>() {
					@Override
					public Boolean call() throws Exception
					{
						try{
							int j = 0;
							double[] rate_array = ga2.crossover_rate(Cross_pool);
							double[] aa = ga2.crossover_rate(ga2.pop);
							
							while(j < ga2.geneSize / 4){
								ga2.crossover(ga2.pop, aa);
								ga2.crossover(Cross_pool, rate_array);
								j++;
							}
							j=0;
							while(j < ga2.geneSize/2){	
								ga2.mutate_SA(ga2.geneSize + j);
								ga2.mutate_SA(ga2.geneSize + j + 1);
								j += 2;
							}
							ga2.select();
							return true;
						}
						catch(Exception e){
							e.printStackTrace();
							return false;
						}
					}
				};
				Future<Boolean> r1;
				Future<Boolean> r2;
				Boolean if_complete1;
				Boolean if_complete2;
				
				ga1.pop.clear();
				ga1.Step = 0;
				ga2.pop.clear();
				ga2.Step = 0;
				Cross_pool.clear();
				
				//thread-	
				Future<Boolean> c1 = es.submit(new Callable<Boolean>() {
					@Override
					public Boolean call() throws Exception
					{
						ga1.createBeginningSpecies();
						return true;
					}
				});
				Future<Boolean> c2 = es.submit(new Callable<Boolean>() {
					@Override
					public Boolean call() throws Exception
					{
						ga2.createBeginningSpecies();
						return true;
					}
				});
				c1.get();
				c2.get();
				
				//put individuals into Cross_pool 
				tools.copy_list(ga1.pop, Cross_pool);
				tools.copy_list(ga2.pop, Cross_pool);
				
				//maxt
				int maxt = 0;
				double pre_fitness = ga1.pop.get(0).fitness[0] > ga2.pop.get(0).fitness[0] ? ga1.pop.get(0).fitness[0] : ga2.pop.get(0).fitness[0];
				double next_fitness = 0.0;
				while(ga1.Step < ga1.genestep)
				{	
					if(maxt == 10){
						break;
					}
					r1 = es.submit(callable1);
					r2 = es.submit(callable2);
					if_complete1 = (Boolean) r1.get();
					if_complete2 = (Boolean) r2.get();
					if(if_complete1 && if_complete2){
						int popcount = ga1.geneSize - 1;
						for(int j = 0; j < 1; j++){
							SpeciesIndividual tempIndividual = tools.copy_SpeciesIndividual(ga1.pop.get(j));
							SpeciesIndividual tempIndividual2 = tools.copy_SpeciesIndividual(ga2.pop.get(j));
							if(tempIndividual.fitness[0] > ga2.pop.get(popcount).fitness[0]){
								ga2.pop.set(popcount, tempIndividual);
								
							}
							if(tempIndividual2.fitness[0] > ga1.pop.get(popcount).fitness[0]){
								ga1.pop.set(popcount, tempIndividual2);
							}
							popcount--;
						}
					}
					
					tools.sort_pop_list(ga1.pop);
					tools.sort_pop_list(ga2.pop);
					tools.copy_list(ga1.pop, Cross_pool);
					tools.copy_list(ga2.pop, Cross_pool);
					update_Cross_pool();
					
					next_fitness = ga1.pop.get(0).fitness[0] > ga2.pop.get(0).fitness[0] ? ga1.pop.get(0).fitness[0] : ga2.pop.get(0).fitness[0];
					if(pre_fitness < next_fitness){
						pre_fitness = next_fitness;
						maxt = 0;
					}
					else{
						maxt++;
					}
					ga1.Step++;
				}
				es.shutdown();
				SpeciesIndividual temp = new SpeciesIndividual();
				temp = tools.copy_SpeciesIndividual(ga1.pop.get(0).fitness[0] > ga2.pop.get(0).fitness[0] ? ga1.pop.get(0):ga2.pop.get(0));
				SpeciesIndividual_index tempmaxgene = new SpeciesIndividual_index(temp,i+1);
				ga1.maxgene.add(tempmaxgene);	

				end = System.currentTimeMillis();
				one_time = Double.valueOf(end - start)/1000;
				time_count += one_time;
				BigDecimal bd = new BigDecimal(one_time).setScale(4, RoundingMode.UP);
				System.out.println("total time：" + bd.doubleValue() + "s\n");
			}	
			SpeciesIndividual_index max = new SpeciesIndividual_index();
			max = tools.searchMax(ga1.maxgene);
			
			//pring results
			BigDecimal bd = new BigDecimal(time_count/step).setScale(4, RoundingMode.UP);
			System.out.print("NO." + max.index + " time, " + "The optimal gene set is：{ ");
			for(int i =  0; i < max.speciesIndividual.chromosome.length; i++)
			{
				System.out.print(ga1.name[max.speciesIndividual.chromosome[i]] + ", ");
			}
			System.out.println("}");
			System.out.println("Fitness: " + max.speciesIndividual.fitness[0]);
			System.out.println("CO(M): " + max.speciesIndividual.fitness[1]);
			System.out.println("ME(M): " + max.speciesIndividual.fitness[2]);
			System.out.println("\n"+step+" times, the average execution time is：" + bd.doubleValue() + "s");
			
			//P value test
			double max_fitness = max.speciesIndividual.fitness[0];
			int correct = 0;
			for(int i = 0; i < P_step; i++){
				int[] chromosome = new int[k];
				int index = ga1.random.nextInt(65535) % (ga1.name_index.length);
				chromosome[0] = index;
				int j = 1;
				for(; j < k;){
					index = ga1.random.nextInt(65535) % (ga1.name_index.length);
					int m = 0;
					for(; m < k; m++){
						if(chromosome[m] == index){
							break;
						}
					}
					if(m == k){
						chromosome[j] = index;
						j++;
					}
				}
				Object[] args = {chromosome, ga1.data_Array.get(0)};
				double[] fitness_result = (double[]) ga1.method.invoke(ga1, args);
				double temp_fitness = fitness_result[0];
				if(max_fitness > temp_fitness){
					correct++;
				}
			}//for 1000
			System.out.println("p-value is: " + (double)correct / P_step);
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}
	
	
	public static void main(String[] args)
	{
		Run r = new Run();
	
		String path = "GBM_removeGene_GeneNumbers_911.txt;";
		String[] paths = path.split(";");
		
		int g = 911;
		int k = 6;
		int size = g / 2;
		
		//The first   parameter: File path，
		//The second  parameter: Number of genes，
		//The third   parameter: Gene set size (k)，
		//The fourth  parameter: Population size (N)，
		//The fifth   parameter: Iteration steps (maxg)，
		//The sixth   parameter: Mutation probability (Pm)，
		//The seventh parameter: Number of times the algorithm is executed，
		//The eighth  parameter: Number of cycles when calculating p-value,
		//The ninth   parameter: Model name ("calfitness_Cov", "calfitness_01"),
		//"calfitness_Cov": Model based on coefficient of variation,
		//"calfitness_01":  The original maximum weight sub-matrix solution model
		r.run(paths, g, k, size, 1000, 0.3, 1, 1000, "calfitness_Cov");
	}
}


