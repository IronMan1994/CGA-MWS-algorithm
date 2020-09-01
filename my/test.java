package my;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.Reader;
import java.lang.reflect.Method;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

public class test
{
	//初始化数据
	public List<Sample> A;
	//基因序列长度
	int n;
	//每个染色体内的基因数
	int k;	
	//种群大小
	int geneSize;
	//每次遗传算法迭代的次数
	int genestep;
	//执行遗传算法的总次数
	int MaxStep;
	//每一步
	int Step;
	//变异概率
	double pm;
	//基因姓名
	public String[] name;
	//基因初始顺序
	public int[] name_index;
	//迭代总数后，放入一次执行过后的遗传算法的最优值
	public List<SpeciesIndividual_index> maxgene;
	//初始种群保存每个染色体的适应度
	public List<SpeciesIndividual> tempPOP;
	//初始种群类
	public List<SpeciesIndividual> pop;
	//下代种群
	public List<SpeciesIndividual> nextPOP;
	//随机变量
	public Random random;
	//计算适应度函数
	public Method method;
	//时间
	public long start;
	public long end;
	
	public String[] path;
	
	//初始化3组数据
	public void initData(String[] paths, int n, int k, int geneSize, int genestep, double pm, int step, String methodName) throws Exception
	{
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
		this.tempPOP = new ArrayList<SpeciesIndividual>();
		this.nextPOP = new ArrayList<SpeciesIndividual>();
		this.random = new Random(System.currentTimeMillis());
		this.path = paths;	
		
		this.A = new ArrayList<Sample>();
		
		File f1 = new File("");
		for(int i = 0; i < paths.length;i++)
		{
			File f = new File(this.getClass().getResource("/" + paths[i]).getPath());
			f1 = f;
			List<float[]> B = new ArrayList<float[]>();
			readData(f, B);
			Sample a = new Sample(B);
			A.add(a);
		}
		
		//读取基因名字
		Reader reader = new FileReader(f1);
		BufferedReader bfr = new BufferedReader(reader);
		String temp;
		try
		{
			while((temp = bfr.readLine()) != null)
			{
				String[] ppp = temp.split("\t");
				if(ppp.length == this.n + 1)
				{
					for(int i = 0; i < this.n; i++)
					{
//						name[i] = ppp[i + 1].replaceAll("[^0-9a-zA-Z\u4e00-\u9fa5.，,。？“”]+","");
						name[i] = ppp[i + 1];
					}
					break;
				}
			}
			for(int i = 0; i < this.n; i++)
			{
				name_index[i] = i;
			}
		} catch (Exception e)
		{
			// TODO: handle exception
		}
		bfr.close();
		reader.close();
		
		//调用哪个适应度计算函数
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
	
	//读取文件函数
	public void readData(File f, List<float[]> B) throws Exception
	{
		InputStream fis = new FileInputStream(f);
		Reader isr = new InputStreamReader(fis);
		//Reader read = new FileReader(f1);
		BufferedReader bfr = new BufferedReader(isr);
		String tempstr;
		int line = 1;
		try
		{
			while((tempstr = bfr.readLine()) != null)
			{
				//分隔符
				String[] ppp = tempstr.split("\t");
				
				if(ppp.length == this.n + 1)
				{	
					if(line == 2)
					{
						float[] temp = new float[this.n];
						for(int i = 0; i < this.n; i++)
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
	//构造初始种群
	public void createBeginningSpecies() throws Exception
	{
//		int count = 0;
//		List<Integer> index = new ArrayList<Integer>();
//		for(int i = 0; i < name_index.length; i++){
//			index.add(i);
//		}
//		while(count < geneSize)
//		{
//			Collections.shuffle(index);
//			//找到name_index每行对应的基因的适应度最大的基因序列，放入pop中
//			int left = 0;
//			tempPOP.clear();
//			for(int l = 0; l < n / k; l++)
//			{
//				//染色体数组
//				int[] chromosome = new int[k];
//				for(int m = 0; m < k; m++)
//				{
//					chromosome[m] = index.get(left + m);
//				}
//				double[] fitness_result;
//				Object[] args = {chromosome, this.A.get(0)};
//				fitness_result = (double[]) method.invoke(this, args);
//				
//				SpeciesIndividual speciesIndividual = new SpeciesIndividual(chromosome, fitness_result);
//				tempPOP.add(speciesIndividual);
//				left += k;
//			}
//				
//			//存放到pop中
//			double max = tempPOP.get(0).fitness[0];
//			SpeciesIndividual speciesIndividual = new SpeciesIndividual();
//			for(int k = 0; k < tempPOP.size(); k++)
//			{
//				
//				if(tempPOP.get(k).fitness[0] >= max)
//				{
//					max = tempPOP.get(k).fitness[0];
//					speciesIndividual = tempPOP.get(k);
//				}
//			}
//			pop.add(speciesIndividual);
//			count++;
//		}//while	
				
		
		
		pop.clear();
		int count = 0;
		List<Integer> index = new ArrayList<Integer>();
		for(int i = 0; i < name_index.length; i++){
			index.add(i);
		}
		while(count < this.geneSize)
		{
			Collections.shuffle(index);
			//前K个留下
			int[] chromosome = new int[k];
			for(int m = 0; m < k; m++)
			{
				chromosome[m] = index.get(m);
			}
			
			double[] fitness_result;
			Object[] args = {chromosome, this.A.get(0)};
			fitness_result = (double[]) method.invoke(this, args);
			SpeciesIndividual speciesIndividual = new SpeciesIndividual(chromosome, fitness_result);
			pop.add(speciesIndividual);
			count++;
		}//while	
		
	}
	
	//轮盘赌(选择)
	public void select()
	{
		nextPOP.clear();
		sort_pop_list(pop);
		
		//最好的
		for(int i = 0; i < geneSize; i++){
			nextPOP.add(pop.get(i));
		}
		copy_list(nextPOP, pop);
		
	}
	
	//轮盘赌(选择)
//	public void select()
//	{
//		nextPOP.clear();
//		//先把适应度最大的放入下一代列表
//		SpeciesIndividual first = new SpeciesIndividual();
//		sort_pop_list();
//		first = pop.get(0);
//		int fuzhiNum = pop.size() / pop.size();
//		for(int i = 0; i < fuzhiNum; i++)
//		{
//			nextPOP.add(first);
//		}
//
//		//轮盘赌选择
//		double[] rate_array = new double[geneSize - fuzhiNum + 1];
//		//计算pop中的总适应度的和
//		double rate = 0.0;
//		double count = sum_fitness(fuzhiNum, pop.size());
//		
//		//System.out.println(count);
//		//把0放进去，为了判断概率的时候好处理
//		rate_array[0] = 0.0;
//
//		for(int i = 1; i < rate_array.length; i++)
//		{
//			rate += pop.get(fuzhiNum - 1 + i).fitness[0] / count;
//			rate_array[i] = rate;
//		}
//		
//		for(int i = 0; i < geneSize - fuzhiNum; i++)
//		{
//			//随机一个0-1的double类型值
//			double r = Math.random();
//			for(int j = 0; j < rate_array.length - 1; j++)
//			{
//				if(r >= rate_array[j] && r <= rate_array[j + 1])
//				{
//					SpeciesIndividual temp = new SpeciesIndividual();
//					temp = pop.get(j);
//					nextPOP.add(temp);
//				}
//			}
//		}
//		copy_list(nextPOP, pop);	
//	}
	
	public double[] crossover_rate(List<SpeciesIndividual> Cross_pool)
	{
		int geneSize = Cross_pool.size();
		int fuzhiNum = 0;
		//轮盘赌选择
		double[] rate_array = new double[geneSize - fuzhiNum + 1];
		//计算pop中的总适应度的和
		double rate = 0.0;
		double sum = sum_fitness(fuzhiNum, geneSize, Cross_pool);
		
		//把0放进去，为了判断概率的时候好处理
		rate_array[0] = 0.0;
		for(int i = 1; i < rate_array.length; i++)
		{
			rate += Cross_pool.get(fuzhiNum - 1 + i).fitness[0] / sum;
			rate_array[i] = rate;
		}
		return rate_array;
	}
	
	//轮盘赌(选择)
	public void crossover(List<SpeciesIndividual> Cross_pool, double[] rate_array) throws Exception
	{
//		sort_pop_list(Cross_pool);
//		
//		//轮盘赌选择
//		double[] rate_array = new double[this.geneSize + 1];
//		//计算pop中的总适应度的和
//		double rate = 0.0;
//		
//		//System.out.println(count);
//		//把0放进去，为了判断概率的时候好处理
//		for(int i = 0; i < this.geneSize; i++)
//		{
//			rate += 2.0 * (this.geneSize - i) / (geneSize * (geneSize + 1));
//			rate_array[i + 1] = rate;
//		}
		

		//随机一个0-1的double类型值
		double r1 = Math.random();
		int index1 = 0;
		SpeciesIndividual parent1 = null;
		for(int j = 0; j < rate_array.length - 1; j++)
		{
			if(r1 >= rate_array[j] && r1 <= rate_array[j + 1])
			{
				parent1 = Cross_pool.get(j);
				index1 = j;
//				System.out.println(j);
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
		
//		//开始交叉
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
//					temp_set.remove(parent2.chromosome[i]);
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
				for(int a : remainList1){
					System.out.print(a + "\t");
				}
				System.out.println();
				
				for(int a : parent2.chromosome){
					System.out.print(a + "\t");
				}
				System.out.println();
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
		
//		List<Integer> temp_list = new ArrayList<Integer>();	
//		for(Integer a : remainList1){
//			temp_list.add(a);
//		}
//		for(Integer a : remainList2){
//			temp_list.add(a);
//		}
//		Collections.shuffle(temp_list);
//		int j = temp_list.size() / 2;
//		for(int i = 0; i < temp_list.size() / 2; i++)
//		{
//			chromosome1[count] = temp_list.get(i);
//			chromosome2[count] = temp_list.get(j++);
//			count++;
//		}
		
		
//		for(int i = 0; i < remainList1.size(); i++)
//		{
//			boolean wi = random.nextBoolean();
//			if(wi == false){
//				chromosome1[count] = remainList1.get(i);
//				chromosome2[count] = remainList2.get(i);
//			}
//			else{
//				chromosome1[count] = remainList2.get(i);
//				chromosome2[count] = remainList1.get(i);
//			}
//			count++;
//		}
		
		
		
		Object[] args = {chromosome1, this.A.get(0)};
		double[] fitness_result = (double[]) method.invoke(this, args);
		SpeciesIndividual speciesIndividual = new SpeciesIndividual(chromosome1, fitness_result);
		
		
		Object[] args2 = {chromosome2, this.A.get(0)};
		double[] fitness_result2 = (double[]) method.invoke(this, args2);
		SpeciesIndividual speciesIndividual2 = new SpeciesIndividual(chromosome2, fitness_result2);
		
		pop.add(speciesIndividual);
		pop.add(speciesIndividual2);
	}
	
	
	
	
	
	//变异
	public void mutate_SA(int index) throws Exception
	{
		double rate = Math.random();
		if(rate <= this.pm){
			mutate(index);
		}
	}
	
	
//	public void mutate(int Spindex) throws Exception
//	{
//		//列表作好删除
//		List<Integer> temp = new ArrayList<Integer>();
//		int[] chromosome = new int[k];
//		chromosome = pop.get(Spindex).chromosome;
//		Object[] args = {chromosome, this.A.get(0)};
//		double init_fitness = ((double[]) method.invoke(this, args))[0];
//		for(int l = 0; l < chromosome.length; l++)
//		{
//			temp.add(chromosome[l]);
//		}
//		int removeIndex = random.nextInt(65535) % k;
//		temp.remove(removeIndex);
//		
//		int[] houxuan = new int[n - k];
//		int[] index = new int[k - 1];
//		//没有删除的基因放入index数组
//		for(int j = 0; j < temp.size(); j++)
//		{
//			index[j] = temp.get(j);
//		}
//		int count = 0;
//		boolean ifcunzai;
//		for(int j = 0; j < name_index.length; j++)
//		{	
//			ifcunzai = false;
//			for(int k = 0; k < chromosome.length; k++)
//			{
//				if(name_index[j] == chromosome[k])
//				{
//					ifcunzai = true;
//					break;
//				}
//			}
//			if(ifcunzai == false)
//			{
//				houxuan[count] = name_index[j];
//				count++;
//			}			
//		}
//		double maxfitness = 0.0;
//		double tempfitness = 0.0;
//		int maxIndex = 0;
//		
//		double[] max_result = new double[3 * this.A.size() + 1];
//		
//
//		for(int k = 0; k < Math.sqrt(name_index.length); k++)
//		{
//			int m = random.nextInt(houxuan.length);
//			double[] tempfitness_result = new double[3 * this.A.size() + 1];	
//			int[] tempArray;
//			temp.add(houxuan[m]);
//			//列表转成数组
//			tempArray = changeListToArray(temp);
//			
//			Object[] args1 = {tempArray, this.A.get(0)};
//			tempfitness_result = (double[]) method.invoke(this, args1);
//			tempfitness = tempfitness_result[0];
//			if(maxfitness < tempfitness){
//				maxfitness = tempfitness;
//				maxIndex = m;
//				max_result = tempfitness_result;
//				
//			}
//			temp.remove(this.k - 1);
//		}
//		//如果找出的适应度比原来的大，则保留，否则基因组不变
//		if(maxfitness > init_fitness){
//			int[] tempArray = new int[k];
//			int maxgene = houxuan[maxIndex];
//			temp.add(this.k - 1, maxgene);
//			tempArray = changeListToArray(temp);	
//			//赋给pop
//			pop.get(Spindex).chromosome = tempArray;
//			pop.get(Spindex).fitness = max_result;
//		}
//		else {
//			pop.get(Spindex).chromosome = chromosome;
//		}
//	}
	
	//变异
	public void mutate(int Spindex) throws Exception
	{
		//列表作好删除
		int[] chromosome = new int[k];
		chromosome = pop.get(Spindex).chromosome;
		double init_fitness = pop.get(Spindex).fitness[0];
		
		//找出k-1适应值大的基因并保留
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
		
//		double max_fitness = -Double.MAX_VALUE;
//		for(int i = 0; i < chromosome.length; i++){
//			int count = 0;
//			int[] temp_chromosome = new int[this.k - 1];
//			for(int j = 0; j < chromosome.length; j++){
//				if(j != i){
//					temp_chromosome[count++] = chromosome[j];
//				}
//			}
//			Object[] args1 = {temp_chromosome, this.A.get(0)};
//			double[] tempfitness_result = (double[]) method.invoke(this, args1);
//			double tempfitness = tempfitness_result[0];
//			if(tempfitness > max_fitness){
//				for(int k = 0; k < temp_chromosome.length; k++){
//					max_gene[k] = temp_chromosome[k];
//				}
//				max_fitness = tempfitness;
//			}
//		}
		
		
		
		//找出候选基因
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
//			int k = random.nextInt(houxuan.length);
//			int k = random.nextInt(houxuan.size());
//			max_gene[this.k-1] = houxuan[k];
			max_gene[this.k-1] = houxuan.get(step);
			Object[] args1 = {max_gene, this.A.get(0)};
			double[] tempfitness_result = (double[]) method.invoke(this, args1);
			tempfitness = tempfitness_result[0];
			if(maxfitness < tempfitness){
				maxfitness = tempfitness;
//				maxIndex = k;
				maxIndex = step;
				max_result = tempfitness_result;
			}
			step++;
		}
		//如果找出的适应度比原来的大，则保留，否则基因组不变
		if(maxfitness > init_fitness){
			max_gene[this.k - 1] = houxuan.get(maxIndex);
			//赋给pop
			pop.get(Spindex).chromosome = max_gene;
			pop.get(Spindex).fitness = max_result;
		}
	}
	
	
	
	//找总迭代次数后中的最优基因
	public SpeciesIndividual_index searchMax()
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
	
	
	public int[] getBest_Bad_SpeciesIndividual_index()
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
	
	
	//计算单个适应度(变异系数模型)
	public double[] calfitness_Cov(int[] chromosome, Sample B)
	{	
		//变异系数模型
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
			int sum_line_half = 0;
			double max = 0;
			for(int j = 0; j < chromosome.length;j++)
			{
				if(temp[chromosome[j]] > max)
				{
					max = temp[chromosome[j]];
				}	
				if(temp[chromosome[j]] == 0){
					sum_line_half++;
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
//					cov = 0.5;
//					cov = max;
//					cov = Math.sqrt(fangcha_fenzi / (this.k-1));
//					if(sum_line_half >= (chromosome.length - 1)){
//						if(this.k <= 6){
////							cov = Math.sqrt(fangcha_fenzi / (this.k-1));
//							cov = 0.5;
//						}
//					}
				}
				else{
//					cov /= (Math.sqrt(this.k)/1.5);
				}
				red += cov;
			}
			ifcancer += max;	
		}
		
//		if(chromosome.length == 2){
//			red /= chromosome.length;
//		}
//		else {
//			red /= chromosome.length - 1;
//		}
		BigDecimal bd = new BigDecimal(ifcancer).setScale(4,RoundingMode.UP);
		ifcancer = bd.doubleValue();
		
		BigDecimal bd2 = new BigDecimal(red).setScale(4,RoundingMode.UP);
		red = bd2.doubleValue();
		
		fitness = (ifcancer + red);
		
		BigDecimal bd3 = new BigDecimal(fitness).setScale(4,RoundingMode.UP);
		fitness = bd3.doubleValue();
		
		
		BigDecimal bd4 = new BigDecimal(ifcancer).setScale(4,RoundingMode.UP);
		
		result[0] = fitness;
		
//		result[1] = ifcancer;
//		result[2] = red;
		
		result[1] = x;
		result[2] = y;
		result[3] = bd4.doubleValue();
		result[4] = Double.valueOf(red);
		
		return result;	
	}
	
	
	
	//计算单个适应度2x-y
	public double[] calfitness_Mine(int[] chromosome, Sample B)
	{
		// 2x-y模型
		double[] result = new double[3];
		double fitness = 0.0;
		double red = 0.0,ifcancer = 0.0;
		List<float[]> sample_one = B.A;	
		for(int i = 0; i < sample_one.size();i++)
		{
			float[] temp = B.A.get(i);
			float max = 0;
			for(int j = 0; j < chromosome.length;j++)
			{
				if(temp[chromosome[j]] > max)
				{
					max = temp[chromosome[j]];
				}	
				red += temp[chromosome[j]];
			}
			ifcancer += max;
		}
		
		fitness = 2 * ifcancer - red;
		result[0] = fitness;
		result[1] = Double.valueOf(ifcancer);
		result[2] = Double.valueOf(red);
		return result;	
			
		
		//0-1模型
//		double[] result = new double[3];
//		double fitness = 0.0;
//		double red = 0.0,ifcancer = 0.0;
//		List<float[]> sample_one = B.A;	
//		for(int i = 0; i < sample_one.size();i++)
//		{
//			float[] temp = B.A.get(i);
//			int tempifcancer = 0;
//			for(int j = 0; j < chromosome.length;j++)
//			{
//				if(temp[chromosome[j]] == 1)
//				{
//					red++;
//					tempifcancer = 1;
//				}
//			}
//			ifcancer += tempifcancer;		
//		}
//		
//		fitness = 2 * ifcancer - red;
//		result[0] = fitness;
//		result[1] = Double.valueOf(ifcancer);
//		result[2] = Double.valueOf(red);
//		return result;
	}
	
	
	//计算单个适应度2x-y(0-1模型)
	public double[] calfitness_01(int[] chromosome, Sample B)
	{
		//0-1模型
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
	
	
	//计算总适应度
	public double[] calfitness(int[] chromosome) throws Exception
	{
		//计算
		double[] result = new double[4 * this.A.size() + 2];
		double[] result_ifcancer = new double[this.A.size()];
		double[] result_red = new double[this.A.size()];
		double[] result_2if_red  = new double[this.A.size()];
		double result_2if_red_fangcha  = 0.0;
		double fitness = 0.0;
		double result_2if_red_avg = 0.0;
		//计算几个样本的适应度
		for(int i = 0; i < this.A.size(); i++)
		{
			Object[] args = {chromosome, this.A.get(i)};
			double[] fitness_one = (double[]) method.invoke(this, args);
			result[i + 1] = fitness_one[0];
			result_ifcancer[i] = fitness_one[1];
			result_red[i] = fitness_one[2];
			result_2if_red[i] = 2*fitness_one[1]-fitness_one[2];
			result_2if_red_avg += result_2if_red[i];
		}
		
		result_2if_red_avg = result_2if_red_avg / this.A.size();
		//求方差
		for(int i = 0; i < this.A.size(); i++)
		{
			result_2if_red_fangcha += Math.pow(result_2if_red[i]- result_2if_red_avg, 2);
		}
		result_2if_red_fangcha = result_2if_red_fangcha / this.A.size();
		result_2if_red_fangcha = Math.sqrt(result_2if_red_fangcha);
		BigDecimal bd1 = new BigDecimal(result_2if_red_fangcha).setScale(5,RoundingMode.UP);
		result[result.length-1]=bd1.doubleValue();

		
		//把ifcaner放入返回数组，然后放回类里面
		for(int i = 0; i < result_ifcancer.length;i++)
		{
			result[this.A.size()+ 1 + i] = result_ifcancer[i];
			result[2*this.A.size() + 1 + i] = result_red[i];
			result[3*this.A.size() + 1 + i] = result_2if_red[i];
		}
				
		//计算总适应度
		for(int i = 0; i < this.A.size(); i++)
		{
			fitness += result[i + 1];
		}
		
		//3个适应度两两差的绝对值
		double absolute = 0.0;
		absolute = Math.abs(result[1] - result[2]) + Math.abs(result[1] - result[3]) + Math.abs(result[2] - result[3]);

		double left = fitness / this.A.size();
		
		
		//平均数
		double avg = fitness / this.A.size();
		double fangcha = 0.0;
		//求方差
		for(int i = 0; i < this.A.size(); i++)
		{
			fangcha += Math.pow(result[i + 1]- avg, 2);
		}
		double right = Math.pow((fangcha / this.A.size()),1.0/1.5);
		//自己模型
		//返回只求和的结果
		//BigDecimal bd = new BigDecimal(fitness).setScale(5,RoundingMode.UP);
		//返回求了方差的结果
		BigDecimal bd = new BigDecimal(left - right).setScale(6,RoundingMode.UP);
		//BigDecimal bd = new BigDecimal(left - absolute).setScale(6,RoundingMode.UP);

		result[0] = bd.doubleValue();
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
		
	
	//遍历总次数中数据的排序(冒泡)
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
	
	
	//计算pop中剩下19个的总适应度
	public double sum_fitness(int fuzhiNum, int size, List<SpeciesIndividual> pop)
	{
		double count = 0.0;
		for(int i = fuzhiNum; i < geneSize; i++)
		{
			count += pop.get(i).fitness[0];
		}
		return count;
	}
	
	
	//写每代的结果
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
			for(int i =0; i < pop.size(); i++)
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
	
	
	//写找到的最优值
	public void writeFileMax(SpeciesIndividual_index temp) throws Exception
	{
		File f = new File("");
		String pathname = f.getCanonicalPath();
		File ff = new File(pathname + "/src/results/result_MAX.txt");
		OutputStream os = new FileOutputStream(ff, true);
		//BufferedWriter bfw = new BufferedWriter(writer);
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
		//bfw.close();
		os.close();
	}
	
	
	//写总迭代的文件
	public void writeFile(List<SpeciesIndividual_index> temp)throws Exception
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
	
	
	//列表转一位数组
	public int[] changeListToArray(List<Integer> list)
	{
		int[] temp = new int[list.size()];
		for(int i = 0; i < list.size(); i++)
		{
			temp[i] = list.get(i);
		}
		return temp;
	}
	
	
	//打印多个样本数据
	public void printSample(List<Sample> A)
	{
		//打印所有样本数据
		for(int i = 0; i < A.size(); i++)
		{
			Sample a = A.get(i);
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
	
	
	//打印一维数组(double)
	public void printArray(double[] temp)
	{
		for(int i = 0; i < temp.length; i++)
		{
			System.out.print(temp[i] + "\t");
		}
	}
	
	
	//拷贝列表
	public void copy_list(List<SpeciesIndividual> a, List<SpeciesIndividual> b)
	{
		b.clear();
		for(int i = 0; i < a.size(); i++)
		{
			SpeciesIndividual temp = copy_SpeciesIndividual(a.get(i));
			b.add(temp);
		}
	}
	
	//拷贝个体
	public SpeciesIndividual copy_SpeciesIndividual(SpeciesIndividual from)
	{
		int[] chromosome = new int[this.k];
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










