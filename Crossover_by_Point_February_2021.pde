//Date: February 4, 2021
//Branch 2: Crossover by Point. Parents exchange genes before and after a crossover point. 
import java.util.ArrayList;
import java.util.Collections;

//Keep track of crossover and mutation episodes
int track_c, track_m;
int e_genes, m_genes;

//Mutation and crossover rates
int crossover_rate = 100;
int mutation_rate = 1;

//Drawing control. If set to "true" phenotypes are drawn. If set to "false" phenotypes are not drawn. 
boolean drawing_activated = false;

//Controls which reproduction scheme is used. If set to 1, reproduction scheme 1 is used (i.e.parents produce one offspring). If set to 2, reproduction scheme 2 is used (i.e. parents produce two offspring) 
int which_scheme = 1;
int reproduction_control;

//On/off switches for fitness criterion
boolean s_criterion1 = true; 
boolean s_criterion2 = true; 
boolean s_criterion3 = true; 
boolean s_criterion4 = true; 
boolean s_criterion5 = true; 
boolean s_criterion6 = true; 
boolean s_criterion7 = true; 
boolean s_criterion8 = true; 
boolean s_criterion9 = true;

// ArrayList for the genomes of the individuals in the population: 
ArrayList<ArrayList<Integer>> genes_pop = new ArrayList<ArrayList<Integer>>(); //Stores the initial population of genotypes. Before evolution starts, genotypes are copied from ArrayList: g_population to ArrayList: genotypes.   
ArrayList<ArrayList> seats_pop = new ArrayList<ArrayList>(); // Contains ArrayLists, one per individual in the population
ArrayList fitness_pop = new ArrayList();// ArrayList of double values to record the fitness of each individual in the population on every generation. To save a history, output to a csv file on every generation 

int generations_max = 1000; //Number of cycles that the algorithm runs through the population 
int g_length = 100; //Length of genotypes
int population = 100; //Size of population
int growth_max = 50; //Maximum growth for a stem (both upwards or downwards) // this could be regulated through fitness evaluation btw
int stem_grid = 1600; //Number of stem cubes in the grid that comprises the phenotype of an individual 
int numOfClusters = 10, sizeMinClusters = 4, sizeMaxClusters = 10;
int[] clusterSizes = new int[numOfClusters]; // An array of (arbitrary) clusters' sizes. These clusters determine upward growth (4 regions) and downward growth (6 regions) -sort of 'organs' of individuals
int[] clusterAnchors = new int[numOfClusters]; // An array of (random) clusters' initial points or 'anchors' determined at setup time   

int t = 20; //Size of the stems 
float rotX, rotY; //Variables to control camera and rotate view

PrintWriter output; // to save files
int[] rSeeds = {021202011, 121202011, 221202011, 321202011, 421202011 }; // an array of int values to replicate runs
String nameRun;

void setup() {
  int thisSeed = rSeeds[0];
  randomSeed(thisSeed); // use int values to replicate runs
  size(1000, 1000, P3D); //
  nameRun = "r"+ thisSeed + "_g" + generations_max + "_l" + g_length + "_p" + population + "_w" + growth_max + "_s" + stem_grid + "_c" + clusterSizes.length;
  // This file is to save initial population for reference:
  output = createWriter("setup_" + nameRun + ".csv");
  genes_pop = new ArrayList<ArrayList<Integer>>();
  //Create initial gene population of i genotypes of j genes each. Genes consist of integer numbers between 0 and growth_max. 
  for (int i = 0; i < population; i++) {
    ArrayList<Integer> genotype = new ArrayList<Integer>(); //Create one genotype per individual in the population 
    output.print(i +":");
    for (int j = 0; j < g_length; j++) {        
      genotype.add(int(random(0, growth_max))); //Select the value and add gene to genotype
      output.print("\t" + genotype.get(j));
    }
    output.print("\n");
    genes_pop.add(genotype); //Add genotype to ArrayList: gene_population  
  } // i initial gene population is defined 
  output.flush();
  output.close(); // closes file "setup_" + nameRun + ".csv"

  //println("Clusters: ");
  for (int c = 0; c < numOfClusters; c++) { // cluster creation done by defining cluster size and initial point or anchor
    clusterSizes[c] = (int)random(sizeMinClusters, sizeMaxClusters);
    //print("clusterSizes[" + c + "] = " + clusterSizes[c] + "\t");  
    clusterAnchors[c] = (int)random(0, (g_length - (clusterSizes[c])) );
    //print("clusterAnchors[" + c + "] = " + clusterAnchors[c] + "\n");
  } // c

  //// Uncomment this to 'debug' cluster creation:
  //output = createWriter("clusters_" + nameRun + ".csv");
  //for (int i = 0; i < population; i++) {
  //  ArrayList g = genes_pop.get(i);
  //  for (int j = 0; j < numOfClusters; j++) {      
  //    int cs = clusterSizes[j], cp = clusterAnchors[j];
  //    output.print("individual " + i + ", cluster " + j + " of size " + cs + " starting in gene " + cp); 
  //    for (int k = 0; k < cs; k++) {
  //      output.print("\t" + g.get(k+cp));
  //    } // k
  //    output.print("\n");
  //  } //j   
  //  output.print("\n");
  //}// i output file to debug cluster creation
  //output.flush();
  //output.close();

  // initial population is set, now run the evolution up to a number of generations_max:
  for (int xs = 0; xs < generations_max; xs++) {
    println("running generation " + xs);    
    decodePhenotypes(xs); // first read genome to build the phenotype of every individual
    evolve(xs); // then evaluate the fitness of each phenotype
    // if (xs = 0) save all images of phenotypes
    // else if (xs == generations_max-1) save all images of phenotypes
  } // xs  

  output = createWriter("final_" + nameRun + ".csv");
  for (int i = 0; i < genes_pop.size(); i++) {
    ArrayList<Integer> genotype = genes_pop.get(i);
    output.print(i +":");
    for (int j = 0; j < g_length; j++) {
      output.print("\t" + genotype.get(j));
    } // j
    output.print("\n");
  }  // i
  output.flush();
  output.close();
  
  if(drawing_activated == true){  
  show_phenotypes();
  }
  //exit();
  
  //for (int i = 0; i < stem_grid; i++) { //? ToDo: I think this is to save one individual to run draw() but let's instead save the data as files and display all final individuals...
  //  //seat_growth_up.add(seats.get(seats.size()-1).get(i).get(0).get(1)); //"49" is the id of the phenotype (from the last generation) that the program will draw / "i" controls the stem accessed each time / "0" gives access to the data set that controls the positive growth of a stem / "1" gives access to the measure of positive growth of a stem 
  //  //seat_growth_down.add(seats.get(seats.size()-1).get(i).get(1).get(1)); //"49" is the id of the phenotype (from the last generation) that the program will draw / "i" controls the stem accessed each time / "1" gives access to the data set that controls the positive growth of a stem / "1" gives access to the measure of positive growth of a stem
  //} // i stem_grid
} //Closes void setup

void draw() {//Draws a selected seat_phenotype from the last generation 
  background(255, 255, 255);
  translate(width,height/2,-width*2);
  rotateX(rotX);
  rotateY(-rotY);
  
  for(int i = 0; i < stem_grid; i++){//Draw seat_phenotype one stem at a time
      ArrayList draw_stem = (ArrayList)seats_pop.get(0).get(i*2);
      int posX = (int)draw_stem.get(1)*t;
      int posY = 0;
      int posZ = (int)draw_stem.get(2)*t;
      int ug = (int)draw_stem.get(3);
      int dg = (int)draw_stem.get(4);
    
      pushMatrix();
      translate(posX,posY,posZ);
      fill(255,0,0);
      smooth();
      box(t);
      popMatrix(); 
    
      for(int j = 1; j <= ug; j++){//Draw the stem's up growth
        pushMatrix();
        translate(posX,posY-(t*j),posZ);
        fill(0,255,0);
        smooth();
        box(t);
        popMatrix();
      }
    
      for(int k = 1; k <= dg; k++){//Draw the stem's down growth
        pushMatrix();
        translate(posX,posY+(t*k),posZ);
        fill(0,0,255);
        smooth();
        box(t);
        popMatrix();
      }
    }//Close Draw seat
} //Closes draw

void mouseDragged() {
  rotY -=(mouseX - pmouseX) * 0.01;
  rotX -=(mouseY - pmouseY) * 0.01;
}


void decodePhenotypes(int xs) { //
  seats_pop = new ArrayList<ArrayList>(); // Resets ArrayList for the pheontypes of all individuals in the population
  for (int i = 0; i < population; i++) { // for each individual in the population...           
    // for each individual seat (i) its genotype needs to be read to decode the genome and translate it into its phenotype
    // we use the cluster information to select the genes that encode each region of stem cubes
    ArrayList<Integer> regionsGrowth = new ArrayList<Integer>(numOfClusters); // the 10 growth regions for the pheontype (growth from grid: 4 regions upwards, 6 downwards)
    for (int c = 0; c < numOfClusters; c++) { // cycle through all clusters to decode the growth in regions ('organs') of this individual seat
      int cs = clusterSizes[c], cp = clusterAnchors[c]; // size and anchor of each cluster 
      //println("cs = " + cs + "\tcp = " + cp + "\tgenes_pop.get(i).size = " + genes_pop.get(i).size()); 
      ArrayList<Integer> thisCluster = new ArrayList(); // to momentarily store and read the genes in the cluster to translate them into growth values
      int sum_alleles = 0;
      for (int a = 0; a < cs; a++) {
        thisCluster.add(genes_pop.get(i).get(a+cp));
        sum_alleles += genes_pop.get(i).get(a+cp);
      } // a
      Collections.sort(thisCluster); // sorts the cluster values, useful to calculate median and pick small to large values in order
      // using mean or median values is too simple, could instead just use that value rather than complicate life with clusters
      //int mean_allele = sum_allele/cs; //int max_allele = Collections.max(genes_in_cluster); // could use average, max, etc...       
      // gene interaction needs to be more interesting, where a single gene can 'flip' a phenotype significantly, so we use modulo https://processing.org/reference/modulo.html
      int allele; // just keep an eye on the sizeMinClusters, sizeMaxClusters values to avoid getting an IndexOutOfBoundsException here:
      if (sum_alleles%10 == 0) allele = thisCluster.get(thisCluster.size()-1); // uses max value only ocassionally, when the sum of alleles is a multiple of 10 <- 'expensive' from a biological viewpoint
      else if (sum_alleles%3 == 0) allele = thisCluster.get(thisCluster.size()/2); //  more common but still rare, median value when the sum is a multiple of 3 
      else if (sum_alleles%2 == 0) allele = thisCluster.get(2); // more common (about 40% of the time) use the third smallest value // 
      else allele = thisCluster.get(0); // and the 'cheaper' most common version is to pick the smallest value in the cluster

      regionsGrowth.add(allele); // assigns the growth value for each region in the grid (0 to 3 are upwards growth, 4 to 9 are downwards growth)
      // regionsGrowth.set(c, median_allele); // assigns the growth value for each region in the grid (0 to 3 are upwards growth, 4 to 9 are downwards growth)
    } // c

    ArrayList<ArrayList> seat_phenotype = new ArrayList<ArrayList>(); // ArrayList for the phenotype for each individual seat
    int gridWidth = int(sqrt(stem_grid));
    int onefourth = int(sqrt(stem_grid)/4);
    // Uncomment this to 'debug' stem cube update:
    //output = createWriter("generation " + xs + " stems of individual " + i + "_" + nameRun + " grid width " + gridWidth + ".csv");
    //output.println("ID\tx\ty\tstem_upgrowth\tstem_downgrowth");    
    for (int s = 0; s < stem_grid; s++) { // for each of the stem cubes in the grid of this individual: 
      ArrayList stem_phenotype = new ArrayList();
      // check to what upward region it belongs (0 to 3) and to what downward region in belongs (4 to 9)
      // for every stem cube, create an ArrayList with: ID, coord, integer of growth up, integer of growth down, boolean for gaps, String for material, etc
      int x = floor(s%gridWidth); // an x to refer to the position of this stem cube (x,y coord)
      int y = floor(s/gridWidth); // an y to refer to the position of this stem cube (x,y coord)            
      int stem_upgrowth =-1, stem_downgrowth =-1; // initialised to -1 for debugging purposes, to make sure all stem cubes belong to a region 
      // checks stem cube's location to define its upgrowth value as defined by the upward clusters: 
      if ( ((x < onefourth) && (y <= onefourth)) || ((x >= gridWidth-onefourth) && (y <= onefourth)) ) { // 
        stem_upgrowth = regionsGrowth.get(0); // region 0 (A and A') is defined by cluster 0
      } else if ( ((x < onefourth) && (y > onefourth)) || ((x >= gridWidth-onefourth) && (y > onefourth)) ) { // 
        stem_upgrowth = regionsGrowth.get(1); // region 1 (B and B') is defined by cluster 1
      } else if ( ((x >= onefourth) && (x < gridWidth-onefourth) && (y <= onefourth))  ) {
        stem_upgrowth = regionsGrowth.get(2); // region 2 (C) is defined by cluster 2
      } else if ( ((x >= onefourth) && (x < gridWidth-onefourth) && (y > onefourth))  ) {
        stem_upgrowth = regionsGrowth.get(3); // region 3 (D) is defined by cluster 3
      } // end of upward growth regions

      // now checks stem cube's location to define its downgrowth value as defined by the downward clusters:
      if ( ((x < onefourth) && (y < onefourth)) || ((x >= gridWidth-onefourth) && (y < onefourth)) ) {
        stem_downgrowth = regionsGrowth.get(4); // region 4 (AA and AA') is defined by cluster 4
      } else if ( ((x < onefourth) && (y >= onefourth) && (y < onefourth*2)) || ((x >= gridWidth-onefourth) && (y >= onefourth) && (y < onefourth*2)) ) {
        stem_downgrowth = regionsGrowth.get(5); // region 5 (BB and BB') is defined by cluster 5
      } else if ( ((x < onefourth) && (y >= onefourth*2) && (y < onefourth*3)) || ((x >= gridWidth-onefourth) && (y >= onefourth*2) && (y < onefourth*3)) ) {
        stem_downgrowth = regionsGrowth.get(6); // region 6 (CC and CC') is defined by cluster 6
      } else if ( ((x < onefourth) && (y >= onefourth*3)) || ((x >= gridWidth-onefourth) && (y >= onefourth*3))   ) {
        stem_downgrowth = regionsGrowth.get(7); // region 7 (DD and DD') is defined by cluster 7
      } else if ( ((x >= onefourth) && (x < gridWidth-onefourth) && (y < onefourth)) || ((x >= onefourth) && (x < gridWidth-onefourth) && (y >= onefourth*3)) ) {
        stem_downgrowth = regionsGrowth.get(8); // region 8 (EE and EE') is defined by cluster 8
      } else if ( ((x >= onefourth) && (x < gridWidth-onefourth) && (y >= onefourth) && (y < onefourth*2)) ||
        ((x >= onefourth) && (x < gridWidth-onefourth) && (y >= onefourth*2) && (y < onefourth*3))) {
        stem_downgrowth = regionsGrowth.get(9); // region 9 (FF and FF') is defined by cluster 9
      }

      if ((stem_downgrowth <= -1) || (stem_downgrowth <= -1) ) println("There is an error defining stem_downgrowth or stem_downgrowth in stem_grid " + s);
      //output.println(s + "\t" + x + "\t" + y + "\t" + stem_upgrowth + "\t" + stem_downgrowth);
      stem_phenotype.add(s); // #0 is the ID of this stem cube
      stem_phenotype.add(x); // #1 is its x position in the grid
      stem_phenotype.add(y); // #2 is its y position in the grid
      stem_phenotype.add(stem_upgrowth); // #3 is its growth upwards
      stem_phenotype.add(stem_downgrowth); // #4 is its growth upwards      
      //boolean shows = false;
      //stem_phenotype.add(shows); // #5 is a swtich on/off to show/hide this stem cube or a section of its growth
      // here could also add variables to control growth of a stem cube sideways at a certain z value

      // when finished defining all stem cube's properties, add it to the ArrayList seat_phenotype
      seat_phenotype.add(stem_phenotype); // #0 of seat_phenotype is an ArrayList of size  stem_grid with all the info per stem cube 
      seat_phenotype.add(regionsGrowth); // #1 of seat_phenotype is an ArrayList of size numOfClusters with the growth values of each region of this seat       
      //color colour = #FFCC00;
      //seat_phenotype.add(colour); // #2 is the colour of this seat 
      //seat_phenotype.add("topmaterial"); // #3 is a material01 of this seat
      //seat_phenotype.add("legsmaterial"); // #4 is a material02 of this seat, etc...
    } //s
    //output.flush();
    //output.close();
    // add Arraylist seat_phenotype to the ArraYlist seats_pop:
    seats_pop.add(seat_phenotype); //
  }// i loop that creates all phenotypes in the population
} // decodePhenotypes()



void show_phenotypes(){//Produces images of all phenotype_seats from last generation  
  background(255, 255, 255);
  translate(width,height/2,-width*2);
  rotateX(-PI/6);
  rotateY(-PI/6);
  rotateZ(-PI);
 
  for(int u = 0; u < seats_pop.size(); u++){
    background(255,255,255);
    for(int i = 0; i < stem_grid; i++){//Draw seat_phenotype one stem at a time
      ArrayList draw_stem = (ArrayList)seats_pop.get(u).get(i*2);
      int posX = (int)draw_stem.get(1)*t;
      int posY = 0;
      int posZ = (int)draw_stem.get(2)*t;
      int ug = (int)draw_stem.get(3);
      int dg = (int)draw_stem.get(4);
    
      pushMatrix();
      translate(posX,posY,posZ);
      fill(255,0,0);
      smooth();
      box(t);
      popMatrix(); 
    
      for(int j = 1; j <= ug; j++){//Draw the stem's up growth
        pushMatrix();
        translate(posX,posY+(t*j),posZ);
        fill(0,255,0);
        smooth();
        box(t);
        popMatrix();
      }
    
      for(int k = 1; k <= dg; k++){//Draw the stem's down growth
        pushMatrix();
        translate(posX,posY-(t*k),posZ);
        fill(0,0,255);
        smooth();
        box(t);
        popMatrix();
      }
    }//Close Draw seat
    String name = "Seat phenotype_" + str(u) +"_" + nameRun;
    save(name);
  }
}//Ends show_phenotypes


void evolve(int xs) { //
  // First, evaluate the fitness of phenotypes:
  //Create .csv file to document fitness. A file is created every time the algorithm covers a 10% of the total number of cycles. 
  
   if(xs == 0 || xs == ((10*generations_max)/100)-1 || xs == ((20*generations_max)/100)-1 || xs == ((30*generations_max)/100)-1 || xs == ((40*generations_max)/100)-1 || xs == ((50*generations_max)/100)-1 || xs == ((60*generations_max)/100)-1 || xs == ((70*generations_max)/100)-1 || xs == ((80*generations_max)/100)-1 || xs == ((90*generations_max)/100)-1 || xs == ((100*generations_max)/100)-1){
     output = createWriter("fitness_" + xs + "_" + nameRun + ".csv");//Add data to .csv file. 
     output.println("r0\tr1\tr2\tr3\tr4\tr5\tr6\tr7\tr8\tr9\tstability\tsize1\tsize2\tseatability\tstdev1\tstdev2\tbackrest\tarmrests\tharmony\tfitness");
   }
   
   fitness_pop = new ArrayList(); // Reset fitness_pop
   for (int i = 0; i < seats_pop.size(); i++) { 
     // for each individual in population...
     ArrayList a_seat = seats_pop.get(i); 
     double fitness = calculateFitness(a_seat);//Add the fintess of phenotye to ArrayList: f_values.
     fitness_pop.add(fitness);
   }
   //Closes and saves the .csv file that stores fitness data
   output.flush();
   output.close();
  
  // Second, create a roulette to select parents:  
  ArrayList<Integer> roulette = new ArrayList<Integer>(); //
  //output = createWriter("generation " + xs + " roulette_" + nameRun + ".csv");
  for (int i = 0; i < fitness_pop.size(); i++) {    
    for (int j = 0; j < (double)fitness_pop.get(i); j++) {
      roulette.add(i); // ArrayList is a very long list of indices that refer to genotypes in genes_pop, proportional to their fitness scores (mentioned as many times)
      //output.print(i + ", ");
    } //j
  } //i  
  //output.flush();
  //output.close();
  //println(roulette);

  // Third, parent selection. 
  ArrayList<Integer> parents = new ArrayList<Integer>(); // Create ArrayList parents a temporary index of size population*2 used to select parents from genes_pop for mating
  while (parents.size() < population*2) {  // ArrayList:parents is an ArrayList of parents selected proportionally to their fitness determined by their # of mentions in roulette 
    int idp1 = (int)random(0,roulette.size());//Select two random locations in roulette.  
    int idp2 = (int)random(0,roulette.size());
   
    int p1 = roulette.get(idp1);//Retrieve the genotype ids contained at those locations and assign them to parent 1 (p1) and parent 2 (p2).
    int p2 = roulette.get(idp2);
    if(p1 != p2){//Compare the id of parent 1 and parent 2. If their ids differ, add them to ArrayList:parents. This ensures the selection of pairs of parents composed of two different genotypes.
      parents.add(roulette.get(idp1));
      parents.add(roulette.get(idp2));
    }
  } // i Parent selection ends
  
  // Fourth, Crossover via random point selection 
  
  //Create .csv file to document crossover episodes. A file is created every time the algorithm covers a 10% of the total number of cycles. 
  if(xs == 0 || xs == ((10*generations_max)/100)-1 || xs == ((20*generations_max)/100)-1 || xs == ((30*generations_max)/100)-1 || xs == ((40*generations_max)/100)-1 || xs == ((50*generations_max)/100)-1 || xs == ((60*generations_max)/100)-1 || xs == ((70*generations_max)/100)-1 || xs == ((80*generations_max)/100)-1 || xs == ((90*generations_max)/100)-1 || xs == ((100*generations_max)/100)-1){
    output = createWriter("crossover_generation_" + xs + "_" + nameRun + ".csv");
  }
  
  //Regulates which reproduction scheme is implemented
  if(which_scheme == 1){
    reproduction_control = population;
  }else if(which_scheme == 2){
    reproduction_control = population/2;
  }
  
  ArrayList<ArrayList<Integer>> offspring_pop = new ArrayList<ArrayList<Integer>>(); //Create ArrayList: offspring_pop. This ArrayList will temporarily store the new generation of genotypes. 
  for (int i = 0; i < reproduction_control; i++) { // for as many times as population size... // in the future this is the point where the size of population could vary depending on average fitness or other factors 
    int p1 = i*2;//Retrieve pairs of parents at specific locations within ArrayList: parents. 
    int p2 = p1+1;
    ArrayList<Integer> parent1 = new ArrayList<Integer>(genes_pop.get(parents.get(p1))); // the genotype of parent 1
    ArrayList<Integer> parent2 = new ArrayList<Integer>(genes_pop.get(parents.get(p2))); // the genotype of parent 2
    
    //Reproduction scheme 1: each pair of parents produces ONE offspring
    if(which_scheme == 1){
    float dice = random(0,100);//Determines if crossover takes place
    if(dice <= crossover_rate){//Crossover occurs
      int cp = round(random(0,g_length-1));//Select a crossover point
      track_c++; //Update the number of crossover episodes that have taken place 
      
      //For each parent, create subLists containing the genes before and after the crossover point. 
      ArrayList<Integer> chop1 = new ArrayList<Integer>(parent1.subList(0,cp));
      ArrayList<Integer> chop2 = new ArrayList<Integer>(parent1.subList(cp,g_length));
      ArrayList<Integer> chop3 = new ArrayList<Integer>(parent2.subList(0,cp));
      ArrayList<Integer> chop4 = new ArrayList<Integer>(parent2.subList(cp,g_length));
      
      e_genes = e_genes + (g_length - chop1.size());
      
      output.print(i + ":" + "\t" + cp + "\t" + chop1.size() + "\t" + chop2.size());//Add data to .csv file. i = offspring id/ cp = crossover point / chop1.size() = no. of genes before cp / chop2.size() = no. of genes after cp
      output.print("\n");
    
      ArrayList<Integer> o1 = new ArrayList<Integer>();//Create ArrayList for offsrping 1 (o1)
      ArrayList<Integer> o2 = new ArrayList<Integer>();//Create ArrayList for offsrping 2 (o2)
      //Add genes to o1
      o1.addAll(chop1);
      o1.addAll(chop4);
      //Add genes to o2
      o2.addAll(chop3);
      o2.addAll(chop2);
    
      float which_o = random(-1,1);//Randomly select whether o1 or o2 is added to offspring_pop
      if(which_o >= -1 && which_o <0){
        offspring_pop.add(o1);
      }else if(which_o >= 0){
        offspring_pop.add(o2);
      }
    }else if(dice > crossover_rate){//Crossover doesn't occur
      ArrayList<Integer> o1 = new ArrayList<Integer>();//Create ArrayList for offsrping
      float dice2 = random(-1,1);//Randonmly select which parent is passed to next generation
      if(dice2 <= 0){
        o1.addAll(parent1);
      }else if(dice2 > 0){
        o1.addAll(parent2);
      }
      offspring_pop.add(o1);
      output.print(i + ":");//Add data to .csv file. i = offspring id
      output.print("\n");
    }
    }//Reproduction scheme 1 closes. 
    
   //Reproduction scheme 2: each pair of parents produces TWO offspring
    if(which_scheme == 2){
      float dice = random(0,100);//Determines if crossover takes place
    
    if(dice <= crossover_rate){//If crossover occurs...
      int cp = round(random(0,g_length-1));//Select a crossover point
      track_c++;
      
      //For each parent, create subLists containing the genes before and after the crossover point. 
      ArrayList<Integer> chop1 = new ArrayList<Integer>(parent1.subList(0,cp));
      ArrayList<Integer> chop2 = new ArrayList<Integer>(parent1.subList(cp,g_length));
      ArrayList<Integer> chop3 = new ArrayList<Integer>(parent2.subList(0,cp));
      ArrayList<Integer> chop4 = new ArrayList<Integer>(parent2.subList(cp,g_length));
      
      e_genes = e_genes + (g_length - chop1.size());
      
      output.print(i + ":" + "\t" + cp + "\t" + chop1.size() + "\t" + chop2.size());//Add data to .csv file. i = offspring id/ cp = crossover point / chop1.size() = no. of genes before cp / chop2.size() = no. of genes after cp
      output.print("\n");
    
      ArrayList<Integer> o1 = new ArrayList<Integer>();//Create ArrayList for offsrping 1 (o1)
      ArrayList<Integer> o2 = new ArrayList<Integer>();//Create ArrayList for offsrping 2 (o2)
      //Add genes to o1
      o1.addAll(chop1);
      o1.addAll(chop4);
      //Add genes to o2
      o2.addAll(chop3);
      o2.addAll(chop2);
      
      offspring_pop.add(o1);//Add offspring 1 to offspring_pop
      offspring_pop.add(o2);//Add offspring 2 to offspring_pop
      
    }else if(dice > crossover_rate){//Crossover doesn't occur...
      offspring_pop.add(parent1);//Add a copy of parent 1 to offspring_pop
      offspring_pop.add(parent2);//add a copy of parent 2 to offspring_pop
      
      output.print(i + ":");//Add data to .csv file. i = offspring id
      output.print("\n");
    } 
    }//Reproduction option 2 closes
    
    //output.print(i +":");
    //output.print("\n        p1: ");
    //for (int q = 0; q < parent1.size(); q++) output.print("\t" + parent1.get(q));
    //output.print("\n        p2: ");
    //for (int q = 0; q < parent2.size(); q++) output.print("\t" + parent2.get(q));
    //output.print("\nchild " + cp + "(" + cs + "): ");
    //for (int q = 0; q < offspring_pop.get(i).size(); q++) output.print("\t" + offspring_pop.get(i).get(q)); // verifies that the offspring is created as expected
    //output.print("\n");
  }// i Crossover ends
  
  //Closes and saves the .csv file that stores crossover data
  output.print("Crossover episodes: " + "\t" + track_c); // Number of crossover episodes that have taken place at cycle "xs"
  output.print("\n");
  output.print("Exchanged genes: " + "\t" + e_genes); // Amount of genes that have been exchanged at cycle "xs"
  output.print("\n");
  
  output.flush();
  output.close();
  
  //Mutation
  
  //Create .csv file to document mutation episodes. A file is created every time the algorithm covers a 10% of the total number of cycles. 
  if(xs == 0 || xs == ((10*generations_max)/100)-1 || xs == ((20*generations_max)/100)-1 || xs == ((30*generations_max)/100)-1 || xs == ((40*generations_max)/100)-1 || xs == ((50*generations_max)/100)-1 || xs == ((60*generations_max)/100)-1 || xs == ((70*generations_max)/100)-1 || xs == ((80*generations_max)/100)-1 || xs == ((90*generations_max)/100)-1 || xs == ((100*generations_max)/100)-1){
    output = createWriter("mutation_generation_" + xs + "_" + nameRun + ".csv");
  }
  
  for (int i = 0; i < offspring_pop.size(); i++) {//For each genotype in ArrayList offspring_pop...
    float mp = random(0,100);//random (-1,1);//random(0, 10000); //Determine if mutation will take place
   
    if (mp <= mutation_rate) {//If mutation takes place...
      int mutation_type = int(floor(random(1,4)));//Determine the type of mutation
      track_m++; // Update the count of mutation episodes
      m_genes = m_genes + mutation_type; // Update the count of mutated genes
      //println("mutation activated generation " + xs + " - genotype: " + i + " - mutation type: " + mutation_type); 
      
      output.print(i + ":" + "\t" + mutation_type); //Add data to .csv file. i = genotype id / mutation_type = number of genes in genotyype "i" that mutate 
      output.print("\n");
      // Implement mutation
      for(int r = 0; r < mutation_type; r++){ 
        int mutation = round(random(0, (g_length-1))); //Randomly select a gene
        int mutated_value = round(random(0, growth_max)); //Randonmly select a new value for that gene
        offspring_pop.get(i).set(mutation, (mutated_value)); //Assign new value to selected gene
        //println("mutation activated at generation " + xs);
      }//Implement mutation closes  
    }else if(mp > mutation_rate){//If mutation doesn't occur
      output.print(i + ":");//Add data to .csv file. i = genotype id
      output.print("\n");
    }
  }//Mutation ends
  
  //Closes and saves the .csv file that stores mutation data
  output.print("Mutation episodes: " + "\t" + track_m); //Number of mutation episodes that have taken place at cycle "xs"
  output.print("\n");
  output.print("Mutated genes: " + "\t" + m_genes); //Amount of genes that have been exchanged at cycle "xs"
  output.print("\n");
 
  output.flush();
  output.close();

  genes_pop.clear();
  genes_pop.addAll(offspring_pop); // updates the population of genotypes with the offspring
  //for (int i = 0; i < population; i++) {//Copy genotypes from ArrayList:offspring to ArrayList:genotypes so the evolutionary loop can restart
  //genotypes.add(offspring_pop.get(i));
  //}
} //Evolve ends


double calculateFitness(ArrayList a_seat) {
  double f = -1;
  // NOTE: for ALL fitness measurements, the range of possible growth is from 0 to growth_max     
  // we read from the region growth to assess fitness criteria based on growth; later add other fitness criteria using colour, materials, etc.
  ArrayList regionsGrowth = (ArrayList)a_seat.get(1); // check in decodePhenotypes(): #1 of seat_phenotype is an ArrayList of size numOfClusters with the growth values of each region of this individual
  int[] regions_g = new int[regionsGrowth.size()]; // create a copy version in int[] format for math calculations below 
  for (int k = 0; k < regionsGrowth.size(); k++) regions_g[k] = (int)regionsGrowth.get(k);    

  // 1st criterion is stability: first check the difference between region 4 (AA and AA') and region 7 (DD and DD') as a more even growth in these regions would make the seat more stable
  int gap_4_7 = growth_max - abs((int)regionsGrowth.get(4) - (int)regionsGrowth.get(7)); // if gap is closer to 0, higher fitness; if closer to growth_max, lower fitness
  //// then check the diff btw regions 5 (BB and BB') and 6 (CC and CC') as this is also a source of stability
  //int gap_5_6 = growth_max - abs((int)regionsGrowth.get(5) - (int)regionsGrowth.get(6)); // if gap is closer to 0, higher fitness; if closer to growth_max, lower fitness
  // a third source of stability is region 4 (AA and AA') with region 8 (EE and EE'), all of these are different ways of scoring high in stability
  //int gap_4_8 = growth_max - abs((int)regionsGrowth.get(4) - (int)regionsGrowth.get(8)); // if gap is closer to 0, higher fitness; if closer to growth_max, lower fitness
  //// the fourth source of stability is region 7 (DD and DD') with region 8 (EE and EE'), we are leaving region 9 (FF and FF') out of stability but will promote it with other fitness criteria
  //int gap_7_8 = growth_max - abs((int)regionsGrowth.get(7) - (int)regionsGrowth.get(8)); // if gap is closer to 0, higher fitness; if closer to growth_max, lower fitness
  int[] stability_scores = { gap_4_7 }; // gap_4_8,gap_5_6,gap_7_8 // I turned off the others because it gets way too easy for chairs to score high in this criterion from the initial generation with random values
  int stability = max(stability_scores); // the staibility criterion is thus the highest number of the possible ways to score on stability
  float stability_n = norm((float)stability, 0, (float)growth_max);
  if(s_criterion1 == false){
    stability = 0;
    stability_n = 0;
  }

  // 2nd criterion is max size: seats with the highest growth in any region score higher
  int size1 = max(regions_g); // the size criterion is thus the highest number of growth in any region (hence there are numOfClusters ways of scoring for size)
  float size1_n = norm((float)size1, 0, (float)growth_max);
  if(s_criterion2 == false){
    size1 = 0;
    size1_n = 0;
  }

  // 3rd criterion is min size: seats with the lowest growth in any region score higher too
  int size2 = min(regions_g); // the size criterion is thus the highest number of growth in any region (hence there are numOfClusters ways of scoring for size)
  float size2_n = norm((float)size2, 0, (float)growth_max);
  if(s_criterion3 == false){
    size2 = 0;
    size2_n = 0;
  }

  // 4th criterion is seatability: smaller growth in region 3 (D) scores higher, as it allows for better seating
  int seatability = growth_max - (int)regionsGrowth.get(3); // region 3 (D) is the 'cushion area'
  float seatability_n = norm((float)seatability, 0, (float)growth_max);
  if(s_criterion4 == false){
    seatability = 0;
    seatability_n = 0;
  }

  // 5th criterion is variance in upside of the seat: to avoid 'blocks' and promote more interesting shapes, reward differences of growth across regions
  // 6th criterion is variance in downside of the seat: to avoid 'blocks' and promote more interesting shapes, reward differences of growth across regions
  double up_mean = 0.0, up_variance = 0.0, up_stdev = 0.0, down_mean = 0.0, down_variance = 0.0, down_stdev = 0.0;
  for (int k = 0; k < 4; k++) { // remember the first four regions determine upgrowth
    up_mean += regions_g[k];
    output.print(regions_g[k] + "\t");
  }
  up_mean /= 4; // regions_g.length;
  for (int k = 4; k < 10; k++) { // remember the last six regions determine downgrowth
    down_mean += regions_g[k];
    output.print(regions_g[k] + "\t");
  }
  down_mean /= 4; // regions_g.length;
  for (int k = 0; k < 4; k++) up_variance += Math.pow(regions_g[k] - up_mean, 2);
  for (int k = 4; k < 10; k++) down_variance += Math.pow(regions_g[k] - down_mean, 2);
  up_variance /= 4; // regions_g.length;
  down_variance /= 6; // regions_g.length;
  up_stdev = Math.sqrt(up_variance); //
  down_stdev = Math.sqrt(down_variance); //
  up_stdev+= ((0 + growth_max)/2);
  down_stdev+= ((0 + growth_max)/2);
  float stdev1_n = norm((float)up_stdev, 0, (float)growth_max);
  float stdev2_n = norm((float)down_stdev, 0, (float)growth_max);
  if(s_criterion5 == false){
    stdev1_n = 0;
  }
  if(s_criterion6 == false){
    stdev2_n = 0;
  }
  
  // 7th criterion is backrest: higher growth in region 2 (C) scores higher, as it provides a better back support
  int backrest = (int)regionsGrowth.get(2); // region 2 (C) is the 'backrest area'
  float backrest_n = norm((float)backrest, 0, (float)growth_max);
  if(s_criterion7 == false){
    backrest = 0;
    backrest_n = 0;
  }

  // 8th criterion is armrest: some growth but not too much is preferred. Peak is 1/3rd of growth_max
  double armrests = growth_max - abs((growth_max*0.33) - (int)regionsGrowth.get(1)); // region 1 (B) is the 'armest area'
  float armrests_n = norm((float)armrests, 0, (float)growth_max);
  if(s_criterion8 == false){
    armrests = 0;
    armrests_n = 0;
  }
  
  // 9th criterion is harmony between 6 previous criteria
  double harmony = (1 - stability_n) + (1 - size1_n) + (1 - size2_n) + (1 - seatability_n) + (1 - stdev1_n) + (1 - stdev2_n) + (1 - backrest_n) + (1 - armrests_n) ; 
  harmony = norm((float)harmony, 0, (float)8); // normalise by the six previous criteria (each 0 to 1)
  if(s_criterion9 == false){
    harmony = 0;
  }

  f = stability_n + size1_n + size2_n + seatability_n + stdev1_n + stdev2_n + backrest_n + armrests_n + harmony; // overall fitness is given by the linear sum of criteria (later on this can be a weighted sum where weights change over time as priorities, tastes, technology changes)
  output.print(stability_n + "\t" + size1_n + "\t" + size2_n + "\t" + seatability_n + "\t" + stdev1_n + "\t" + stdev2_n + "\t" + backrest_n + "\t" + armrests_n + "\t" + harmony + "\t" + f + "\n" );

  return f;
} // calculateFitness
