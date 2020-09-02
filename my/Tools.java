package my;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.Reader;
import java.util.ArrayList;
import java.util.List;

public class Tools
{
	public Tools()
	{
		
	}
	
	//initData
	public void initData(String[] paths, List<Data_Array> data_Array, int n, String[] name, int[] name_index) throws Exception
	{
		File f = new File("");
		for(int i = 0; i < paths.length;i++)
		{
			File f1 = new File(this.getClass().getResource("/" + paths[i]).getPath());
			f = f1;
			List<float[]> temp_data_array = new ArrayList<float[]>();
			readData(f, temp_data_array, n);
			Data_Array a = new Data_Array(temp_data_array);
			data_Array.add(a);
		}
		
		//读取基因名字
		Reader reader = new FileReader(f);
		BufferedReader bfr = new BufferedReader(reader);
		String temp;
		try
		{
			while((temp = bfr.readLine()) != null)
			{
				String[] ppp = temp.split("\t");
				if(ppp.length == n + 1)
				{
					for(int i = 0; i < n; i++)
					{
						name[i] = ppp[i + 1];
					}
					break;
				}
			}
			for(int i = 0; i < n; i++)
			{
				name_index[i] = i;
			}
		} catch (Exception e)
		{
			e.printStackTrace();
		}
		bfr.close();
		reader.close();
	}
	
	//Read data
	public void readData(File f, List<float[]> B, int n) throws Exception
	{
		InputStream fis = new FileInputStream(f);
		Reader isr = new InputStreamReader(fis);
		BufferedReader bfr = new BufferedReader(isr);
		String tempstr;
		int line = 1;
		try
		{
			while((tempstr = bfr.readLine()) != null)
			{
				//分隔符
				String[] ppp = tempstr.split("\t");
				if(ppp.length == n + 1)
				{	
					if(line == 2)
					{
						float[] temp = new float[n];
						for(int i = 0; i < n; i++)
						{
							//提取数字		
							temp[i] = Float.valueOf(ppp[i + 1].trim());		
						}
						B.add(temp);
					}						
				}
				line = 2;
			}		
		} 
		catch (Exception e)
		{
			e.printStackTrace();
		}
		bfr.close();
		fis.close();
	}
	
	//Find the optimal gene set after the total number of iterations
	public SpeciesIndividual_index searchMax(List<SpeciesIndividual_index> maxgene)
	{
		SpeciesIndividual_index temp = new SpeciesIndividual_index();
		double maxfitness = 0.0;
		double tempfitness = 0.0;
		int maxIndex = 0;
		for(int i = 0; i < maxgene.size(); i++)
		{
			if(maxfitness < maxgene.get(i).speciesIndividual.fitness[0])
			{
				tempfitness = maxgene.get(i).speciesIndividual.fitness[0];
				maxfitness = tempfitness;
				maxIndex = i;
			}
		}
		temp = maxgene.get(maxIndex);
		return temp;
	}
	
	//Sorting of data in total traversal times (bubble algorithm)
	public void sort_list_maxgene(List<SpeciesIndividual_index> temppop)
	{
		SpeciesIndividual_index temp = new SpeciesIndividual_index();
        int size = temppop.size();
        for(int i = 0 ; i < size-1; i ++)
        {
	        for(int j = 0 ;j < size-1-i ; j++)
	        {
	        	//<表示把最小值移到最后
	            if(temppop.get(j).speciesIndividual.fitness[0] < temppop.get(j+1).speciesIndividual.fitness[0])  //交换两数位置
	            {
		            temp = temppop.get(j);
		            temppop.set(j, temppop.get(j+1));
		            temppop.set(j+1, temp);
	            }
	        }
        }
	}
	
	//Population fitness ranking (bubble algorithm)
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
	
	//The function of roulette probability calculation
	public double sum_fitness(int fuzhiNum, int size, List<SpeciesIndividual> pop, int geneSize)
	{
		double count = 0.0;
		for(int i = fuzhiNum; i < geneSize; i++)
		{
			count += pop.get(i).fitness[0];
		}
		return count;
	}
	
	//Find the best and worst individuals of the current population
	public int[] getBest_Bad_SpeciesIndividual_index(List<SpeciesIndividual> pop)
	{
		int[] result = new int[2];
		int max = 0;
		int min = 0;
		
		for(int i = 0; i < pop.size(); i++)
		{
			if(pop.get(max).fitness[0] < pop.get(i).fitness[0]){
				max = i;
			}
			if(pop.get(min).fitness[0] > pop.get(i).fitness[0]){
				min = i;
			}
		}
		
		result[0] = max;
		result[1] = min;
		return result;
	}
	
	//Write the population information for each iteration
	public void writeFileStep(List<SpeciesIndividual> temp, int count) throws Exception
	{
		File f = new File("");
		String pathname = f.getCanonicalPath();
		File ff = new File(pathname + "/src/result_Step.txt");
		if(count == 0)
		{
			ff.delete();
			ff.createNewFile();
		}
		OutputStream os = new FileOutputStream(ff,true);
		try
		{
			os.write(("第" + count + "代、种群信息:").getBytes());
			os.write(("\r\n").getBytes());
			for(int i =0; i < temp.size(); i++)
			{
				for(int j = 0; j < temp.get(i).chromosome.length; j++)
				{
					os.write((String.format("%-8s", temp.get(i).chromosome[j])).getBytes());
				}
				for(int j = 0; j < temp.get(i).fitness.length - 1;j++)
				{
					os.write((String.format("%-12s", temp.get(i).fitness[j])).getBytes());
				}
				os.write("\r\n".getBytes());
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		finally
		{
			os.close();
		}
	}
	
	//Output the optimal gene set found
	public void writeFileMax(SpeciesIndividual_index temp, String[] path, String[] name) throws Exception
	{
		File f = new File("");
		String pathname = f.getCanonicalPath();
		File ff = new File(pathname + "/src/results/result_MAX.txt");
		OutputStream os = new FileOutputStream(ff, true);
		try
		{
			os.write(path[0].getBytes());
			os.write("\r\n".getBytes());
			os.write(("第" + temp.index + "次遗传算法、" + "最优基因组:").getBytes());
			for(int i = 0; i < temp.speciesIndividual.chromosome.length; i++)
			{
				os.write((String.format("%-10s", name[temp.speciesIndividual.chromosome[i]])).getBytes());
			}
			os.write("\r\n".getBytes());
			os.write(("样本分别适应度:").getBytes());
			for(int i = 0; i < temp.speciesIndividual.fitness.length - 1; i++)
			{
				os.write((String.format("%-12s", temp.speciesIndividual.fitness[i + 1])).getBytes());
			}
			os.write("\r\n".getBytes());
			os.write(("最大适应度:" + temp.speciesIndividual.fitness[0]).getBytes());
			os.write("\r\n".getBytes());
			os.write("\r\n".getBytes());
			System.out.println("输出成功!");
		} catch (Exception e){
			e.printStackTrace();
		}
		os.close();
	}

	//Output the total iteration result
	public void writeFile(List<SpeciesIndividual_index> temp, List<SpeciesIndividual_index> maxgene, String[] name)throws Exception
	{
		sort_list_maxgene(maxgene);
		File f = new File("");
		String pathname = f.getCanonicalPath();
		File ff = new File(pathname + "/src/result.txt");
		OutputStream os = new FileOutputStream(ff);
		//OutputStreamWriter os1 = new OutputStreamWriter(os);
		//BufferedWriter bfw = new BufferedWriter(os1);
		try
		{
			for(int i =0; i < maxgene.size(); i++)
			{
				os.write(String.format("第" + "%-2d" + "次遗传算法、" + "最优基因组:",i + 1).getBytes());
				for(int j = 0; j < temp.get(i).speciesIndividual.chromosome.length; j++)
				{
					os.write((String.format("%-10s", name[temp.get(i).speciesIndividual.chromosome[j]])).getBytes());
				}
				os.write("\r\n".getBytes());
				os.write("样本分别适应度:".getBytes());
				for(int j = 0; j < temp.get(i).speciesIndividual.fitness.length - 1;j++)
				{
					os.write((String.format("%-12s", temp.get(i).speciesIndividual.fitness[j + 1])).getBytes());
				}
				os.write("\r\n".getBytes());
				os.write("样本总适应度:".getBytes());
				os.write((temp.get(i).speciesIndividual.fitness[0] +"").getBytes());
				os.write("\r\n".getBytes());
				os.write("\r\n".getBytes());
			}
			
			System.out.println("输出成功!");
		} catch (Exception e){
			e.printStackTrace();
		}
		//bfw.close();
		os.close();
	}
	
	//Record operation log
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
		
	//Deep copy list
	public void copy_list(List<SpeciesIndividual> a, List<SpeciesIndividual> b)
	{
		b.clear();
		for(int i = 0; i < a.size(); i++)
		{
			SpeciesIndividual temp = copy_SpeciesIndividual(a.get(i));
			b.add(temp);
		}
	}
	
	//Copy and add individuals
	public void copy_add_list(List<SpeciesIndividual> a, List<SpeciesIndividual> b)
	{
		for(int i = 0; i < a.size(); i++)
		{
			SpeciesIndividual temp = copy_SpeciesIndividual(a.get(i));
			b.add(temp);
		}
	}
	
	//Deep copy individual
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

	//Convert the list into a one-dimensional array
	public int[] changeListToArray(List<Integer> list)
	{
		int[] temp = new int[list.size()];
		for(int i = 0; i < list.size(); i++)
		{
			temp[i] = list.get(i);
		}
		return temp;
	}
		
	//Print input data matrix
	public void printSample(List<Data_Array> A)
	{
		//打印所有样本数据
		for(int i = 0; i < A.size(); i++)
		{
			Data_Array a = A.get(i);
			for(int j = 0; j < a.A.size(); j++)
			{
				float[] temp = a.A.get(j);
				for(int k1 = 0; k1 < temp.length; k1++)
				{
					System.out.print(temp[k1] + "\t");
				}
				System.out.println();
			}
			System.out.println("\n\n\n");
		}
	}
	
	//Print one dimensional array(double)
	public void printArray(double[] temp)
	{
		for(int i = 0; i < temp.length; i++)
		{
			System.out.print(temp[i] + "\t");
		}
	}

	//Print population information quickly
	public void printPOP(List<SpeciesIndividual> tempPOP)
	{
		for(int i = 0; i < tempPOP.size(); i++)
		{
			for(int j = 0; j < tempPOP.get(i).chromosome.length; j++)
			{			
				System.out.print(String.format("%-12s", tempPOP.get(i).chromosome[j]));
			}
			for(int j = 0; j < tempPOP.get(i).fitness.length; j++)
			{
				System.out.print(String.format("%-12s", tempPOP.get(i).fitness[j]));
			}
			System.out.println();
		}
	}
}
