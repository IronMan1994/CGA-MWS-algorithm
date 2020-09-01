package my;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
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
	//时间
	public long start;
	public long end;
	public List<SpeciesIndividual> Cross_pool = new ArrayList<SpeciesIndividual>();
	
	
	//记录
	public void writeResult(String file_name, List temp, String methond, int k) throws IOException
	{
		File f = new File("");
		String pathname = f.getCanonicalPath();
		File ff = new File(pathname + "/src/results/Results.txt");
		OutputStream os = new FileOutputStream(ff, true);
		try{
			os.write((file_name + "    methond_name = " + methond + "    K = " + k + "：").getBytes());
			os.write("\r\n".getBytes());
			for(int i = 0; i < temp.size(); i++){
				os.write(String.format("%-12s",  " ").getBytes());
				os.write((String.valueOf(temp.get(i))).getBytes());
				os.write("\r\n".getBytes());
			}
			os.write("\r\n\r\n\r\n".getBytes());
		}
		catch (Exception e){
			e.printStackTrace();
		}
		finally{
			os.close();
		}
	}
	
	
	//拷贝列表
	public void copy_list(List<SpeciesIndividual> a, List<SpeciesIndividual> b)
	{
//		b.clear();
		for(int i = 0; i < a.size(); i++)
		{
			SpeciesIndividual temp = copy_SpeciesIndividual(a.get(i));
			b.add(temp);
		}
	}
	
	public void copy_list2(List<SpeciesIndividual> from, List<SpeciesIndividual> to)
	{
		to.clear();
		for(int i = 0; i < from.size(); i++)
		{
			SpeciesIndividual temp = copy_SpeciesIndividual(from.get(i));
			to.add(temp);
		}
	}
	
	public SpeciesIndividual copy_SpeciesIndividual(SpeciesIndividual from)
	{
		int[] chromosome = new int[from.chromosome.length];
		double[] fitness = new double[from.fitness.length];
		
		for(int i = 0; i < chromosome.length; i++){
			chromosome[i] = from.chromosome[i];
		}
		
		for(int i = 0; i < fitness.length; i++){
			fitness[i] = from.fitness[i];
		}
		SpeciesIndividual result = new SpeciesIndividual(chromosome, fitness);
		
		return result;
	}
	
	//适应度排序(冒泡)
	public void sort_pop_list(List<SpeciesIndividual> pop)
	{
		SpeciesIndividual temp = new SpeciesIndividual();
        int size = pop.size();
        for(int i = 0 ; i < size-1; i ++)
        {
	        for(int j = 0 ;j < size-1-i ; j++)
	        {
	        	//<表示把最小值移到最后
	            if(pop.get(j).fitness[0] < pop.get(j+1).fitness[0])  //交换两数位置
	            {
		            temp = pop.get(j);
		            pop.set(j, pop.get(j+1));
		            pop.set(j+1, temp);
	            }
	        }
	        
        }
	}
	
	public void update_Cross_pool(test t)
	{
		List<SpeciesIndividual> tempList = new ArrayList<SpeciesIndividual>();
//		double[] rate_array = t.crossover_rate(Cross_pool);
//		for(int i = 0; i < Cross_pool.size() / 2; i++){
//			double r1 = Math.random();
//			for(int j = 0; j < rate_array.length - 1; j++)
//			{
//				if(r1 >= rate_array[j] && r1 <= rate_array[j + 1])
//				{
//					SpeciesIndividual temp = t.copy_SpeciesIndividual(Cross_pool.get(j));
//					tempList.add(temp);
//				}
//			}
//		}
		
		
		sort_pop_list(Cross_pool);
		for(int i = 0; i < Cross_pool.size() / 2; i++){
			SpeciesIndividual temp = copy_SpeciesIndividual(Cross_pool.get(i));
			tempList.add(temp);
		}
		copy_list2(tempList, Cross_pool);
	}

	
	public void run(String[] paths, int n ,int k, int geneSize, int genestep, double pm,int step, int P_step, String methodName)
	{
		test t = new test();
		test t2 = new test();
	
		try
		{
			//读文件，和初始化全局变量
			t.initData(paths, n, k, geneSize/2, genestep, pm, step, methodName);
			t2.initData(paths, n, k, geneSize/2, genestep, pm, step, methodName);
			//迭代
			List results = new ArrayList<>();
			double one_time = 0.0;
			double time_count = 0.0;
			for(int i = 0; i < t.MaxStep; i++)
			{
				start = System.currentTimeMillis();	
				System.out.println(String.format("第%-2d次执行遗传算法，重新构造初始种群", i + 1));
				ExecutorService es = Executors.newFixedThreadPool(2);
				Callable<Boolean> callable1 = new Callable<Boolean>() {
					@Override
					public Boolean call() throws Exception
					{
						try{
							int j = 0;
							double[] rate_array = t.crossover_rate(Cross_pool);
							double[] aa = t.crossover_rate(t.pop);
							
							while(j < t.geneSize / 4){
								t.crossover(t.pop, aa);
								t.crossover(Cross_pool, rate_array);
								j++;
							}
							j=0;
							while(j < t.geneSize/2)
							{	
//								t.crossover(Cross_pool, rate_array);
								t.mutate_SA(t.geneSize + j);
								t.mutate_SA(t.geneSize + j + 1);
								j += 2;
							}
							t.select();
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
//							while(j < t2.geneSize)
//							{	
//								t2.crossover(t2.geneSize + j);
//								t2.mutate_SA(t2.geneSize + j, 1);
//								j += 1;
//							}
							double[] rate_array = t2.crossover_rate(Cross_pool);
							double[] aa = t2.crossover_rate(t2.pop);
							
							while(j < t2.geneSize / 4){
								t2.crossover(t2.pop, aa);
								t2.crossover(Cross_pool, rate_array);
								j++;
							}
							j=0;
							while(j < t2.geneSize/2)
							{	
//								t2.crossover(Cross_pool, rate_array);
								t2.mutate_SA(t2.geneSize + j);
								t2.mutate_SA(t2.geneSize + j + 1);
								j += 2;
							}
							t2.select();
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
				
				t.pop.clear();
				t.Step = 0;
				t2.pop.clear();
				t2.Step = 0;
				Cross_pool.clear();
//				t.maxgene.clear();
				
				//创建初始种群	
				Future<Boolean> c1 = es.submit(new Callable<Boolean>() {
					@Override
					public Boolean call() throws Exception
					{
						t.createBeginningSpecies();
						return true;
					}
				});
				
				Future<Boolean> c2 = es.submit(new Callable<Boolean>() {
					@Override
					public Boolean call() throws Exception
					{
						t2.createBeginningSpecies();
						return true;
					}
				});
				c1.get();
				c2.get();
				
				//初代种群放入交叉池
				copy_list(t.pop, Cross_pool);
				copy_list(t2.pop, Cross_pool);
				
				//适应度10次不变则跳出
				int count = 0;
				//取出上次迭代最好适应度
				double pre_fitness = t.pop.get(0).fitness[0] > t2.pop.get(0).fitness[0] ? t.pop.get(0).fitness[0] : t2.pop.get(0).fitness[0];
//				double pre_fitness = t.pop.get(0).fitness[0];
				double next_fitness = 0.0;
				while(t.Step < t.genestep)
				{	
					if(count == 10){
						break;
					}
					
					r1 = es.submit(callable1);
					r2 = es.submit(callable2);
					
					if_complete1 = (Boolean) r1.get();
					if_complete2 = (Boolean) r2.get();
					if(if_complete1 && if_complete2){
						int popcount = t.geneSize - 1;
						for(int j = 0; j < 1; j++){
							SpeciesIndividual tempIndividual = t.copy_SpeciesIndividual(t.pop.get(j));
							SpeciesIndividual tempIndividual2 = t2.copy_SpeciesIndividual(t2.pop.get(j));
							
							if(tempIndividual.fitness[0] > t2.pop.get(popcount).fitness[0]){
								t2.pop.set(popcount, tempIndividual);
								
							}
							if(tempIndividual2.fitness[0] > t.pop.get(popcount).fitness[0]){
								t.pop.set(popcount, tempIndividual2);
							}
							popcount--;
						}
					}
					
					t.sort_pop_list(t.pop);
					t2.sort_pop_list(t2.pop);
					
					copy_list(t.pop, Cross_pool);
					copy_list(t2.pop, Cross_pool);
					
					if(t.Step > 1 && t.Step%2 == 0){
						update_Cross_pool(t);
					}
//					update_Cross_pool(t);
					
					//十步不变就跳出
					next_fitness = t.pop.get(0).fitness[0] > t2.pop.get(0).fitness[0] ? t.pop.get(0).fitness[0] : t2.pop.get(0).fitness[0];
//					next_fitness = t.pop.get(0).fitness[0];
//					if(pre_fitness == next_fitness){
//						count++;
//					}
					if(pre_fitness < next_fitness){
						pre_fitness = next_fitness;
						count = 0;
					}
					else{
						count++;
					}
					t.Step++;
				}
				es.shutdown();
				SpeciesIndividual temp = new SpeciesIndividual();
				temp = copy_SpeciesIndividual(t.pop.get(0).fitness[0] > t2.pop.get(0).fitness[0] ? t.pop.get(0):t2.pop.get(0));
//				temp = t.pop.get(0);
				SpeciesIndividual_index tempmaxgene = new SpeciesIndividual_index(temp,i+1);
				t.maxgene.add(tempmaxgene);	
				
				end = System.currentTimeMillis();
				one_time = Double.valueOf(end - start)/1000;
				time_count += one_time;
				BigDecimal bd = new BigDecimal(one_time).setScale(4, RoundingMode.UP);
				
				String str = "";
				str += String.format("%-20s", (i+1) + "：" + String.format("%.04f", bd.doubleValue()) + "秒");;
				for(int j =  0; j < temp.chromosome.length; j++){
					str += String.format("%-10s", t.name[temp.chromosome[j]]);
				}
				str += temp.fitness[0] + " ";
				results.add(str);
				
				System.out.println("执行总时间为：" + bd.doubleValue() + "秒\n");
//				for(int j = 0; j < temp.chromosome.length; j++){
//					System.out.print(temp.chromosome[j] + "\t");
//				}
//				System.out.println();
			}
//			end = System.currentTimeMillis();
			//写总迭代结果
			t.writeFile(t.maxgene);		
			SpeciesIndividual_index max = new SpeciesIndividual_index();
			max = t.searchMax();
			//输出文件
			t.writeFileMax(max);
			//打印结果
			BigDecimal bd = new BigDecimal(time_count/step).setScale(4, RoundingMode.UP);
			System.out.println("第" + max.index + "次遗传算法、" + "最优基因：");
			for(int i =  0; i < max.speciesIndividual.chromosome.length; i++)
			{
				System.out.print(t.name[max.speciesIndividual.chromosome[i]] + "\t");
			}
			System.out.println(max.speciesIndividual.fitness[0] + "\t" + max.speciesIndividual.fitness[1] + "\t" + max.speciesIndividual.fitness[2]);	
			System.out.println("分别适应度为:");
			for(int i = 0; i < max.speciesIndividual.fitness.length - 1; i++)
			{
				System.out.print(max.speciesIndividual.fitness[i + 1] + "\t");
			}
			
			String str = (step + 1)+ "：" + bd.doubleValue() + "秒\t";;
			results.add(str);
			writeResult(paths[0], results, methodName, k);
			
			System.out.println("\n"+step+"次执行平均时间为：" + bd.doubleValue() + "秒");
			
			for(int j = 0; j < t.maxgene.size(); j++){
				for(int i =  0; i < t.maxgene.get(j).speciesIndividual.chromosome.length; i++)
				{
					System.out.print(t.name[t.maxgene.get(j).speciesIndividual.chromosome[i]] + "\t");
				}
				System.out.println();
			}
			
			System.out.println(t.Step);
			
			//测P值
//			double max_fitness = max.speciesIndividual.fitness[0];
//			int correct = 0;
//			
//			for(int i = 0; i < P_step; i++){
//				//随机算k个基因测适应值
//				int[] chromosome = new int[t.k];
//				int index = t.random.nextInt(65535) % (t.name_index.length);
//				chromosome[0] = index;
//				int j = 1;
//				for(; j < t.k;){
//					index = t.random.nextInt(65535) % (t.name_index.length);
//					int m = 0;
//					for(; m < t.k; m++){
//						if(chromosome[m] == index){
//							break;
//						}
//					}
//					//表示没有重复
//					if(m == t.k){
//						chromosome[j] = index;
//						j++;
//					}
//				}
//				double temp_fitness = t.calfitness(chromosome)[0];
//				if(max_fitness > temp_fitness){
//					correct++;
//				}
//			}//for 1000
//			System.out.println("P值为：" + (double)correct / P_step);

		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}
	
	
	public static void main(String[] args)
	{
		Run r = new Run();
		
//		String[] paths = {"I1", "I5", "I10", "G1000I1","G1000I5","G1000I10","G5000I1","G5000I5","G5000I10","G10000I1","G10000I5","G10000I10"};
//		for(int i = 0; i < paths.length; i++){
//			File f = new File("src//"+paths[i]);
//			String files[] = f.list();
//			for(String ss: files)
//			{
//				String path = paths[i] + "//" + ss + ";";
//				String[] pathss = path.split(";");
//				int g = Integer.valueOf(ss.substring(ss.indexOf("G")+1, ss.lastIndexOf("I")-1));
//				int k = 2;
////				int size = (int) (Math.log(Math.pow(g, k)) / Math.log(2)) * 10;
//				int size = g / 2; 
//				System.out.println(path);
//				//前1个是文件路径，
//				//第2个参数是基因数，
//				//第3个参数是，一个染色体中的基因个数，
//				//第4个参数是种群大小，   
//				//第5个是每次遗传算法迭代的次数
//				//第6个是变异概率，
//				//第7个是总的执行遗传算法的次数，
//				//第8个是计算P时候的循环次数
//				//第9个是计算适应值函数(calfitness_Cov，calfitness_Mine, calfitness_01)
//				r.run(pathss, g, k, size, 500, 0.3, 10, 1000, "calfitness_01");
//			}
//		}
//		
//		
//		
//		
//		String path = "G10000I10//simulate_P1000_G10000_I10_K10.txt;";
//		String[] pathss = path.split(";");
//		int g = Integer.valueOf(path.substring(path.lastIndexOf("G")+1, path.lastIndexOf("I")-1));
//		
//		for(int k = 2; k <=10; k++){
////			int size = (int) (Math.log(Math.pow(g, k)) / Math.log(2))*10;
//			int size = g / 2;
//			//前1个是文件路径，
//			//第2个参数是基因数，
//			//第3个参数是，一个染色体中的基因个数，
//			//第4个参数是种群大小，   
//			//第5个是每次遗传算法迭代的次数
//			//第6个是变异概率，
//			//第7个是总的执行遗传算法的次数，
//			//第8个是计算P时候的循环次数
//			//第9个是计算适应值函数(calfitness_Cov，calfitness_Mine, calfitness_01)
//			r.run(pathss, g, k, size, 500, 0.3, 10, 1000, "calfitness_01");
//		}
//		
//		
//		String pathsss = "G10000I10//simulate_P1000_G10000_I10_K10.txt;";
//		String[] pathssss = pathsss.split(";");
//		
//		for(int k = 2; k <=10; k++){
////			int size = (int) (Math.log(Math.pow(g, k)) / Math.log(2))*10;
//			int size = g / 2;
//			//前1个是文件路径，
//			//第2个参数是基因数，
//			//第3个参数是，一个染色体中的基因个数，
//			//第4个参数是种群大小，   
//			//第5个是每次遗传算法迭代的次数
//			//第6个是变异概率，
//			//第7个是总的执行遗传算法的次数，
//			//第8个是计算P时候的循环次数
//			//第9个是计算适应值函数(calfitness_Cov，calfitness_Mine, calfitness_01)
//			r.run(pathssss, g, k, size, 500, 0.3, 10, 1000, "calfitness_Cov");
//		}
		
		
		
////		String path = "G10000I10/simulate_P1000_G10000_I10_K10_0.5.txt;";
//		String path = "A_cut_off-1-sample_cut_off-0_GBM.txt;";
//		String[] paths = path.split(";");
//		
////		int g = Integer.valueOf(path.substring(path.lastIndexOf("G")+1, path.lastIndexOf("I")-1));
//		int g = 920;
//		int k = 10;
//		
////		int size = (int) (Math.log(Math.pow(g, k)) / Math.log(2)) * 10 * 2;
//		
////		int size = (int) Math.sqrt(g*k*30);
//		
////		int size = (int) (g * (k/10.0));
//		
//		int size = g / 2;
//		
//		System.out.println(size);
////		int size = (int) (Math.sqrt(g));
////		int size = (int) Math.sqrt(g)*2;
//		//前1个是文件路径，
//		//第2个参数是基因数，
//		//第3个参数是，一个染色体中的基因个数，
//		//第4个参数是种群大小，   
//		//第5个是每次遗传算法迭代的次数
//		//第6个是变异概率，
//		//第7个是总的执行遗传算法的次数，
//		//第8个是计算P时候的循环次数
//		//第9个是计算适应值函数(calfitness_Cov，calfitness_Mine, calfitness_01)
//		r.run(paths, g, k, size, 500, 0.3, 10, 1000, "calfitness_Cov");
		
		
		
		String[] paths = {"OV_test"};
		for(int i = 0; i < paths.length; i++){
			File f = new File("src//"+paths[i]);
			String files[] = f.list();
			for(String ss: files)
			{
				String path = paths[i] + "//" + ss + ";";
				String[] pathss = path.split(";");
				int g = Integer.valueOf(ss.substring(ss.lastIndexOf("e")+2, ss.lastIndexOf("e")+6));
				int k = 4;
//				int size = (int) (Math.log(Math.pow(g, k)) / Math.log(2)) * 10;
				int size = g / 2; 
				System.out.println(path);
				System.out.println(g);
				
				for(;k <= 6; k++)
				{
					//前1个是文件路径，
					//第2个参数是基因数，
					//第3个参数是，一个染色体中的基因个数，
					//第4个参数是种群大小，   
					//第5个是每次遗传算法迭代的次数
					//第6个是变异概率，
					//第7个是总的执行遗传算法的次数，
					//第8个是计算P时候的循环次数
					//第9个是计算适应值函数(calfitness_Cov，calfitness_Mine, calfitness_01)
					r.run(pathss, g, k, size, 500, 0.3, 10, 1000, "calfitness_Cov");
				}				
			}
		}	
	}
}


