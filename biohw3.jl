include("flux.jl")


using LinearAlgebra


S_matrix = Array{Float64,2}(undef,18,26)
S_matrix[1,:]=[-1	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];#aspar
S_matrix[2,:]=[1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];#arginsucc
S_matrix[3,:]=[0	1	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];#fumarate
S_matrix[4,:]=[0	1	-1	0	-1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];#arginine
S_matrix[5,:]=[0	0	1	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];#urea
S_matrix[6,:]=[0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];#ornithine
S_matrix[7,:]=[-1	0	0	1	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];#citruline
S_matrix[8,:]=[0	0	0	-1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];#carbyphos
S_matrix[9,:]=[-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0];#ATP
S_matrix[10,:]=[1	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0];#diphos
S_matrix[11,:]=[1	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];#AMP
S_matrix[12,:]=[0	0	0	1	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0];#phosphate
S_matrix[13,:]=[0	0	0	0	-1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1]; #oxygen
S_matrix[14,:]=[0	0	0	0	-1	1	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0]; #NADPH
S_matrix[15,:]=[0	0	0	0	-1	1	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0]; #hydrogen
S_matrix[16,:]=[0	0	0	0	1	-1	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	1	0	0	0	0]; #NO
S_matrix[17,:]=[0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	1	0	0	0]; #NADP
S_matrix[18,:]=[0	0	0	0	1	-1	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	1	0	0	0	0	0]; #water




# Kcat values for v's
k_cat=Array{Float64,2}(undef,6,1)
k_cat[1]=203;#v1
k_cat[2]=34.5;#v2
k_cat[3]=249;#v3
k_cat[4]=88.1;#v4
k_cat[5]=13.7;#v5a
k_cat[6]=13.7;#v5b
E=0.01;
theta=1;

#values from paper (or assigned value from Ecoli or 0)  because they werent there for the 18 metabolites
metabolite_value=Array{Float64,2}(undef,18,1)
metabolite_value[1]=0.01492*1000000;#aspar
metabolite_value[2]=0.0;#argisuc
metabolite_value[3]=0.00048*1000000;#fumarate
metabolite_value[4]=0.00025*1000000;#arginine
metabolite_value[5]=0.0;#urea
metabolite_value[6]=0.00448*1000000;#ornithine
metabolite_value[7]=0.00059*1000000#carbyphos - Ecoli carbamoyl-aspartate
metabolite_value[8]=0.0;#citruline
metabolite_value[9]= 0.00006*1000000#ATP
metabolite_value[10]=0.0;#diphos
metabolite_value[11]=0.0000423016355975951*1000000;#AMP
metabolite_value[12]=0.0000284*1000000;#phosphate - ribose-5-phosphate
metabolite_value[13]=0.0;#oxygen
metabolite_value[14]=0.00006502026004935134*1000000;#NADPH
metabolite_value[15]=0.0;#hydrogen
metabolite_value[16]=0.0#NO
metabolite_value[17]=0.000502026004935134*1000000;#NADP
metabolite_value[18]=0.0;#water

cellvolume=1e-12;
cellmass=2.3e-9;
percentwater=0.9;
cellmassdry=(1-percentwater)*cellmass;
convertfactor(value)=value*(cellvolume)/cellmassdry;
metabolite=convertfactor(metabolite_value)

#Km values from paper or made up
kM_value=Array{Float64,2}(undef,11,1)
kM_value[1]=0.000154276895439494*1000000;#aspar
kM_value[2]=0.0000*1000000;#arginsucc
kM_value[3]=0.0053154608643830104*1000000;#fumarate
kM_value[4]=3.49694159840394;#arginine
kM_value[5]=0.0053*1000000;#carbyphos
kM_value[6]=0.00085*1000000;#ornithine
kM_value[7]=0.000056*1000000;#citruline
kM_value[8]=0.0;#citruline
kM_value[9]=0.000392333154234199*1000000;#ATP
kM_value[10]=0.00003*1000000;#NADPH
kM_value[11]=0.000027*1000000#NADP
kM=convertfactor(kM_value)



Default_bounds_array = zeros(26,2)
#Upper bounds
Default_bounds_array[1,2]=k_cat[1]*E*theta*(metabolite[1]/(metabolite[1]+kM[1]))*(metabolite[9]/(metabolite[9]+kM[9]))*1;# Some terms assumed to be at saturation
Default_bounds_array[2,2]=k_cat[2]*E*theta*1;# assumed 1 km = 0
Default_bounds_array[3,2]=k_cat[3]*E*theta*(metabolite[4]/(metabolite[4]+kM[4]));
Default_bounds_array[4,2]=k_cat[4]*E*theta*(metabolite[6]/(metabolite[6]+kM[6]))*(metabolite[7]/(metabolite[7]+kM[5]));#Some terms assumed to be at saturation
Default_bounds_array[5,2]=k_cat[5]*E*theta*(metabolite[4]/(metabolite[4]+kM[4]))*(metabolite[14]/(metabolite[14]+kM[10]));#Some terms assumed to be at saturation
Default_bounds_array[6,2]=k_cat[6]*E*theta*(metabolite[17]/(metabolite[17]+kM[11]));
Default_bounds_array[7,2]=(2.77);
Default_bounds_array[8,2]=(2.77);
Default_bounds_array[9,2]=(2.77);
Default_bounds_array[10,2]=(2.77);
Default_bounds_array[11,2]=(2.77);
Default_bounds_array[12,2]=(2.77);
Default_bounds_array[13,2]=(2.77);
Default_bounds_array[14,2]=(2.77);
Default_bounds_array[15,2]=(2.77);
Default_bounds_array[16,2]=(2.77);
Default_bounds_array[17,2]=(2.77);
Default_bounds_array[18,2]=(2.77);
Default_bounds_array[19,2]=(2.77);
Default_bounds_array[20,2]=(2.77);
Default_bounds_array[21,2]=(2.77);
Default_bounds_array[22,2]=(2.77);
Default_bounds_array[23,2]=(2.77);
Default_bounds_array[24,2]=(2.77);
Default_bounds_array[25,2]=(2.77);
Default_bounds_array[26,2]=(2.77);

#Set lower bounds on auxilliary metabolites as these reactions are reversible
Default_bounds_array[11,1]=-(2.77);
Default_bounds_array[12,1]=-(2.77);
Default_bounds_array[13,1]=-(2.77);
Default_bounds_array[14,1]=-(2.77);
Default_bounds_array[15,1]=-(2.77);
Default_bounds_array[16,1]=-(2.77);
Default_bounds_array[17,1]=-(2.77);
Default_bounds_array[18,1]=-(2.77);
Default_bounds_array[19,1]=-(2.77);
Default_bounds_array[20,1]=-(2.77);
Default_bounds_array[21,1]=-(2.77);
Default_bounds_array[22,1]=-(2.77);
Default_bounds_array[23,1]=-(2.77);
Default_bounds_array[24,1]=-(2.77);
Default_bounds_array[25,1]=-(2.77);
Default_bounds_array[26,1]=-(2.77);


#Create no bounds for 18 metabolites
Species_bounds_array=zeros(18,2);


Objective_coefficient_array=zeros(20);
Objective_coefficient_array[10]=-1;

Min_flag = true;

#Use the FBA function
answer=calculate_optimal_flux_distribution(S_matrix,Default_bounds_array,Species_bounds_array,Objective_coefficient_array);
#Extract the flux vector
r=3600*answer[2]/1000; #mmol/gDW-hr

urea_flux=r[10];
println("MaxFlux",urea_flux)
