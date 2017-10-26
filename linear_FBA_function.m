function [Flux_model,bound_LFBA,solution_FBA,solution_LFBA]=linear_FBA_function()
%% read model and cofactor

cofactor={'adp[c]','atp[c]','amp[c]','nadp[c]','nadph[c]','nad[c]','nadh[c]','h2o[c]','h2o[p]','h[c]','h[e]','h[p]','coa[c]','ACP[c]','pi[c]','pi[p]','ppi[c]','cmp[c]','q8h2[c]','q8[c]','mqn8[c]','mql8[c]','gtp[c]','gdp[c]','ctp[c]','pi[e]','udp[c]','fadh2[c]','gmp[c]','ump[c]','fe2[c]','co2[c]','co2[e]','co2[p]','o2[c]'};
cofactor_hash = java.util.Hashtable;
for i=1:length(cofactor)
    cofactor_hash.put(cofactor{i},i);
end
Ecoli_Model = readCbModel('.\data\Ec_iAF1260_flux1.xml');
%% reaction constraint
exchange_reaction_index=[];
input_nutrient={};
index=1;
secrated_metabolite_set={'SUCCtex','ETOHtex','FORtex'};
% secrated_metabolite_set={''};
secrated_metabolite_positive={'ACtex'};
% secrated_metabolite_positive={''};
input_nutrient_file=fopen('.\data\uptake_nutrient_glucose.txt');
while(~feof(input_nutrient_file))
    line=fgetl(input_nutrient_file);
    if(~isempty(line))
        input_nutrient{index}=line;
        index=index+1;
    end
end
fclose(input_nutrient_file);
input_output_nutrient={};
index=1;
input_nutrient_file=fopen('.\data\general_input_output.txt');
while(~feof(input_nutrient_file))
    line=fgetl(input_nutrient_file);
    if(~isempty(line))
        input_output_nutrient{index}=line;
        index=index+1;
    end
end
fclose(input_nutrient_file);

input_amino_acid={};
index=1;
input_nutrient_file=fopen('.\data\amino_input.txt');
while(~feof(input_nutrient_file))
    line=fgetl(input_nutrient_file);
    if(~isempty(line))
        input_amino_acid{index}=line;
        index=index+1;
    end
end
fclose(input_nutrient_file);

input_amino_acid_limit={};
index=1;
input_nutrient_file=fopen('.\data\amino_input_limited.txt');
while(~feof(input_nutrient_file))
    line=fgetl(input_nutrient_file);
    if(~isempty(line))
        input_amino_acid_limit{index}=line;
        index=index+1;
    end
end
fclose(input_nutrient_file);

for i=1:length(Ecoli_Model.rxns)
    if ~(isempty(regexp(Ecoli_Model.rxns{i},'\w+tex$')) && isempty(regexp(Ecoli_Model.rxns{i},'\w+texi$')))
%         disp(Ecoli_Model.rxns{i})
        flag=0;
        exchange_reaction_index=[exchange_reaction_index,i];
        for j=1:length(secrated_metabolite_set)
            if strcmp(Ecoli_Model.rxns{i},secrated_metabolite_set{j})
                Ecoli_Model.lb(i)=-5;
                Ecoli_Model.ub(i)=2;
                flag=1;
                break;
            end
        end
        if flag==1
            continue;
        end
        for j=1:length(secrated_metabolite_positive)
            if strcmp(Ecoli_Model.rxns{i},secrated_metabolite_positive{j})
%                 Ecoli_Model.lb(i)=-640-10;
%                 Ecoli_Model.ub(i)=-640+10;
                Ecoli_Model.ub(i)=-5;
                flag=1;
                break;
            end
        end
        if flag==1
            continue;
        end
        for j=1:length(input_nutrient)
            if strcmp(Ecoli_Model.rxns{i},input_nutrient{j})
%                 Ecoli_Model.rev(i)=0;
                Ecoli_Model.lb(i)=-1;
%                   Ecoli_Model.lb(i)=5;
                flag=1;
                break;
            end
        end
        if flag==1
            continue;
        end
        
        for j=1:length(input_amino_acid)
            if strcmp(Ecoli_Model.rxns{i},input_amino_acid{j})
%                 Ecoli_Model.rev(i)=0;
                Ecoli_Model.lb(i)=-1;
                Ecoli_Model.ub(i)=50;
%                   Ecoli_Model.lb(i)=5;
                flag=1;
                break;
            end
        end
        if flag==1
            continue;
        end
        for j=1:length(input_amino_acid_limit)
            if strcmp(Ecoli_Model.rxns{i},input_amino_acid_limit{j})
%                 Ecoli_Model.rev(i)=0;
                Ecoli_Model.lb(i)=-1;
                Ecoli_Model.ub(i)=10;
%                   Ecoli_Model.lb(i)=5;
                flag=1;
                break;
            end
        end
        if flag==1
            continue;
        end
        
        for j=1:length(input_output_nutrient)
            if strcmp(Ecoli_Model.rxns{i},input_output_nutrient{j})
%                 Ecoli_Model.lb(i)=5;
                flag=1;
                break;
            end
        end
        if flag==1
            continue;
        end
%         Ecoli_Model.rev(i)=0;
        Ecoli_Model.ub(i)=1;
    end
end
%% construct condition specific network based on nutrient
lowFlux_reaction=[];
initCobraToolbox();
for i=1:length(Ecoli_Model.rxns)
     Ecoli_Model.c(:)=0;
     Ecoli_Model.c(i)=1;
     solution1=optimizeCbModel(Ecoli_Model);
     Ecoli_Model.c(i)=-1;
     solution2=optimizeCbModel(Ecoli_Model);
     if solution1.f==0 && solution2.f==0
         lowFlux_reaction=[lowFlux_reaction,i];
     end
end
Ecoli_Model.c(1:length(Ecoli_Model.c))=0;
specific_reaction=1:length(Ecoli_Model.rxns);
specific_reaction=setdiff(specific_reaction,lowFlux_reaction);
for i=1:length(specific_reaction)
    Ecoli_Model_specific.rxns{i}=Ecoli_Model.rxns{specific_reaction(i)};
end
Ecoli_Model_specific.mets=Ecoli_Model.mets;
Ecoli_Model_specific.S=Ecoli_Model.S(:,specific_reaction);
Ecoli_Model_specific.rev=Ecoli_Model.rev(specific_reaction);
Ecoli_Model_specific.lb=Ecoli_Model.lb(specific_reaction);
Ecoli_Model_specific.ub=Ecoli_Model.ub(specific_reaction);
Ecoli_Model_specific.c=Ecoli_Model.c(specific_reaction);
Ecoli_Model_specific.metCharge=Ecoli_Model.metCharge;
for i=1:length(specific_reaction)
    Ecoli_Model_specific.rules{i}=Ecoli_Model.rules{specific_reaction(i)};
end
Ecoli_Model_specific.genes=Ecoli_Model.genes;
for i=1:length(specific_reaction)
    Ecoli_Model_specific.grRules{i}=Ecoli_Model.grRules{specific_reaction(i)};
end
for i=1:length(specific_reaction)
    Ecoli_Model_specific.subSystems{i}=Ecoli_Model.subSystems{specific_reaction(i)};
end
for i=1:length(specific_reaction)
    Ecoli_Model_specific.rxnNames{i}=Ecoli_Model.rxnNames{specific_reaction(i)};
end
Ecoli_Model_specific.metNames=Ecoli_Model.metNames;
Ecoli_Model_specific.metFormulas=Ecoli_Model.metFormulas;
Ecoli_Model_specific.b=Ecoli_Model.b;
%% construct connect graph
network=ConnectGraphConstruct(Ecoli_Model_specific.S,Ecoli_Model_specific.rxns,Ecoli_Model_specific.mets,cofactor_hash,Ecoli_Model_specific.rev);
%% calculate the simple cycle
flag_java=simplecycle_java();
similarity_matrix=load('similarity_matrix_5len_rec3.txt');
% similarity_matrix=load('similarity_matrix_reac7Length.txt');
%% cluster analysis
 nonzeroReactionSet=[];
 noClusterReactionSet=[];
 for i=1:size(similarity_matrix,1)
     if ~(sum(similarity_matrix(:,i))==0 && sum(similarity_matrix(i,:))==0)
         nonzeroReactionSet=[nonzeroReactionSet,i];
     else
         noClusterReactionSet=[noClusterReactionSet,i];
     end
 end
 
  similarityMatrix_nozero=zeros(length(nonzeroReactionSet));
 for i=1:length(nonzeroReactionSet)
     for j=1:length(nonzeroReactionSet)
         similarityMatrix_nozero(i,j)=similarity_matrix(nonzeroReactionSet(i),nonzeroReactionSet(j));
     end
 end

  cluster_num=5;
  clusterStru=cluster_similarityMatrix(similarityMatrix_nozero,cluster_num);
  cluster_set={};
  for i=1:cluster_num
      temp_index=find(clusterStru==i);
      cluster_set{i}=nonzeroReactionSet(temp_index);
  end
  %% construct sub metabolic model for each sub cluster
  Stichimetic_Submatrix_set={};
  reaction_set=1:size(Ecoli_Model_specific.S,2);
  for i=1:length(cluster_set)
      temp=Ecoli_Model_specific.S(:,cluster_set{i});
      zero_rowIndex=[];
      rest_reaction_set=setdiff(reaction_set,cluster_set{i});
      for j=1:size(temp,1)
          if sum(abs(temp(j,:)))==0
              zero_rowIndex=[zero_rowIndex,j];
          end
      end
      metabolite_inset=1:size(temp,1);
      metabolite_set=setdiff(metabolite_inset,zero_rowIndex);
      temp=temp(metabolite_set,:);
      
      reversibleReac_index=find(Ecoli_Model_specific.rev(cluster_set{i})==1);
      temp1=cluster_set{i};
      temp1=temp1(reversibleReac_index);
      temp2=Ecoli_Model_specific.S(metabolite_set,temp1)*-1;
      temp=[temp2';temp'];
      temp=temp';
      
      
      addtional_colmon=[];
      index=1;
      for j=1:length(metabolite_set)
          input_temp=find(Ecoli_Model_specific.S(metabolite_set(j),rest_reaction_set) > 0);
          input_test=find(temp(j,:) < 0);
          if ~(isempty(input_temp) || isempty(input_test))
            temp_colmon=zeros(1,length(metabolite_set));
            temp_colmon(j)=1;
            addtional_colmon=[addtional_colmon;temp_colmon];
          end
          
          output_temp=find(Ecoli_Model_specific.S(metabolite_set(j),rest_reaction_set) < 0);
          output_test=find(temp(j,:) > 0);
          if ~(isempty(output_temp) || isempty(output_test))
            temp_colmon=zeros(1,length(metabolite_set));
            temp_colmon(j)=-1;
            addtional_colmon=[addtional_colmon;temp_colmon];
          end
      end
      temp=[temp';addtional_colmon];
      
%       reversibleReac_index=find(Ecoli_Model_specific.rev(cluster_set{i})==1);
%       temp1=cluster_set{i};
%       temp1=temp1(reversibleReac_index);
%       temp2=Ecoli_Model_specific.S(metabolite_set,temp1)*-1;
%       temp=[temp2';temp];
%       temp=temp';
      Stichimetic_Submatrix_set{i,1}=temp';
      Stichimetic_Submatrix_set{i,2}=metabolite_set;
      Stichimetic_Submatrix_set{i,3}=cluster_set{i};
      Stichimetic_Submatrix_set{i,4}=temp1;
      reversible_pair_index=[];
      reac_index=cluster_set{i};
      for j=1:length(reversibleReac_index)
          reversible_pair_index(j,1)=j;
          reversible_pair_index(j,2)=length(reversibleReac_index)+reversibleReac_index(j);
      end
      Stichimetic_Submatrix_set{i,5}=reversible_pair_index;
  end
  %% solving the sparse basis vector of sub cluster model based on the L1 approprate method

  for i=1:size(Stichimetic_Submatrix_set,1)
      sub_model(i).numRxns=size(Stichimetic_Submatrix_set{i,1},2);
      sub_model(i).obj=zeros(sub_model(i).numRxns,1);
      sub_model(i).S=sparse(Stichimetic_Submatrix_set{i,1});
      sub_model(i).A=sparse(Stichimetic_Submatrix_set{i,1});
      sub_model(i).rhs=zeros(size(sub_model(i).S,1),1);
      sub_model(i).InternalReacNum= length(Stichimetic_Submatrix_set{i,3})+ size(Stichimetic_Submatrix_set{i,5},1);
      sub_model(i).IOReacNum=sub_model(i).numRxns-sub_model(i).InternalReacNum;
      
      sub_model(i).lb(1:length(Stichimetic_Submatrix_set{i,4}))=0;
      sub_model(i).lb(1+length(Stichimetic_Submatrix_set{i,4}):length(Stichimetic_Submatrix_set{i,3})+length(Stichimetic_Submatrix_set{i,4}))=Ecoli_Model_specific.lb(Stichimetic_Submatrix_set{i,3});
      sub_model(i).lb(length(sub_model(i).lb)+1:sub_model(i).numRxns)=0;
      sub_model(i).lb(find(sub_model(i).lb<0))=0;
      sub_model(i).lb=sub_model(i).lb';
      
      sub_model(i).ub(1:length(Stichimetic_Submatrix_set{i,4}))=abs(Ecoli_Model_specific.lb(Stichimetic_Submatrix_set{i,4}));
      sub_model(i).ub(1+length(Stichimetic_Submatrix_set{i,4}):length(Stichimetic_Submatrix_set{i,3})+length(Stichimetic_Submatrix_set{i,4}))=Ecoli_Model_specific.ub(Stichimetic_Submatrix_set{i,3});
      sub_model(i).ub(length(sub_model(i).ub)+1:sub_model(i).numRxns)=1000;
      sub_model(i).ub=sub_model(i).ub';
      
      temp_pair_reac=Stichimetic_Submatrix_set{i,5};
      temp_binaryMatrix=zeros(size(temp_pair_reac,1),sub_model(i).numRxns);
      for j=1:size(temp_pair_reac,1)
          temp_binaryMatrix(j,temp_pair_reac(j,1))=1;
          temp_binaryMatrix(j,temp_pair_reac(j,2))=1;
      end
      sub_model(i).binaryMatrix=temp_binaryMatrix;
      
      
      sub_model(i).vtype='C';
      sub_model(i).modelsense='max';
      sub_model(i).sense='=';
      
      temp_index=Stichimetic_Submatrix_set{i,4};
      for j=1:length(temp_index)
          sub_model(i).rxns{j}=strcat('rev_',Ecoli_Model_specific.rxns{temp_index(j)});
          sub_model(i).rxnNames{j}=strcat('rev_',Ecoli_Model_specific.rxnNames{temp_index(j)});
      end
      temp_index=Stichimetic_Submatrix_set{i,3};
      for j=1:length(temp_index)
          sub_model(i).rxns{j+length(Stichimetic_Submatrix_set{i,4})}=Ecoli_Model_specific.rxns{temp_index(j)};
          sub_model(i).rxnNames{j+length(Stichimetic_Submatrix_set{i,4})}=Ecoli_Model_specific.rxnNames{temp_index(j)};
      end
      temp_index_meta=Stichimetic_Submatrix_set{i,2};
      for j=length(sub_model(i).rxns)+1:sub_model(i).numRxns
          metaName=find(sub_model(i).S(:,j)~=0);
          if length(metaName)==1
              if sub_model(i).S(metaName(1),j)>0
                  sub_model(i).rxns{j}=strcat('Input_',Ecoli_Model_specific.mets{temp_index_meta(metaName(1))});
                  sub_model(i).rxnNames{j}=strcat('Input_',Ecoli_Model_specific.metNames{temp_index_meta(metaName(1))});
              else 
                  sub_model(i).rxns{j}=strcat('Output_',Ecoli_Model_specific.mets{temp_index_meta(metaName(1))});
                  sub_model(i).rxnNames{j}=strcat('Output_',Ecoli_Model_specific.metNames{temp_index_meta(metaName(1))});
              end
          end
      end
      temp_index=Stichimetic_Submatrix_set{i,2};
      for j=1:length(temp_index)
          sub_model(i).mets{j}=Ecoli_Model_specific.mets{temp_index(j)};
          sub_model(i).metNames{j}=Ecoli_Model_specific.metNames{temp_index(j)};
      end
      sub_model(i).description=strcat('subnetwork_',num2str(i));
      sub_model(i).rev=zeros(sub_model(i).numRxns,1);

      sub_model(i).c(1:length(Stichimetic_Submatrix_set{i,4}))=0;
      sub_model(i).c(1+length(Stichimetic_Submatrix_set{i,4}):length(Stichimetic_Submatrix_set{i,3})+length(Stichimetic_Submatrix_set{i,4}))=Ecoli_Model.c(Stichimetic_Submatrix_set{i,3});
      for j=length(sub_model(i).c)+1:sub_model(i).numRxns
          sub_model(i).c(j)=0;
      end
  end
   % solving the sparse absis vector using the SNP techenique
  Nsnp={};
  for i=1:length(sub_model)
    Nsnp{i} = fastSNP(sub_model(i),'gurobi');
  end
  %% reconstruct new stoichimetic matrix for the sparse basis vector and the rest linear reactions in orginal network
  reaction_set_basisVector={};
  reaction_set_basisVector_lb = [];
  reaction_set_basisVector_ub = [];
  basisVector_set=[];
  reaction_basisVector_index={};
  reaction_basisVector_coff={};
  repeat_hash=java.util.Hashtable;
  index=1;
  for i=1:length(Nsnp)
      current_basis_vectors=Nsnp{i};
      
      metabolite_index=Stichimetic_Submatrix_set{i,2};
      temp_reactionset=[];
      internal_reaction_num=sub_model(i).InternalReacNum;
      
      reaction_set_index=[];
      temp_rev=Stichimetic_Submatrix_set{i,5};
      temp_rev=temp_rev(:,2);
      temp_rev=temp_rev-length(temp_rev);
      temp_rev_1=Stichimetic_Submatrix_set{i,3};
      temp_rev=temp_rev_1(temp_rev);
      reaction_set_index=[temp_rev,temp_rev_1];
      reaction_set_rev=zeros(1,length(reaction_set_index));
      reaction_set_rev(1:length(temp_rev))=1;
      
      for j=1:size(current_basis_vectors,2)
          temp_basisVector=current_basis_vectors(:,j);
          temp_basisVector=temp_basisVector(1:internal_reaction_num);
          temp_metabolites=sub_model(i).S(:,1:internal_reaction_num)*temp_basisVector;
          if sum(temp_metabolites(find(temp_metabolites ~= 0)))<1.0e-8
              single_reaction=find(temp_basisVector ~= 0);
              temp_reaction=[];
              for k=1:length(single_reaction)
                  if repeat_hash.isEmpty || ~repeat_hash.containsKey(reaction_set_index(single_reaction(k)))
                      if single_reaction(k)<=length(temp_rev)
                          repeat_hash.put(reaction_set_index(single_reaction(k)),-1);
                          temp_reaction=[temp_reaction,single_reaction(k)];
                      else
                          repeat_hash.put(reaction_set_index(single_reaction(k)),1);
                          temp_reaction=[temp_reaction,single_reaction(k)];
                      end
                  else
                      if single_reaction(k)<=length(temp_rev) && repeat_hash.get(reaction_set_index(single_reaction(k)))==1
                          temp_reaction=[temp_reaction,single_reaction(k)];
                          repeat_hash.put(reaction_set_index(single_reaction(k)),0);
                      elseif single_reaction(k)> length(temp_rev) && repeat_hash.get(reaction_set_index(single_reaction(k)))==-1
                          temp_reaction=[temp_reaction,single_reaction(k)];
                          repeat_hash.put(reaction_set_index(single_reaction(k)),0);
                      end
                  end
              end
              if isempty(temp_reaction)
                  continue;
              else
                  single_reaction=temp_reaction;
              end
              temp_lb=sub_model(i).lb(single_reaction);
              temp_ub=sub_model(i).ub(single_reaction);
              reaction_set_basisVector_lb=[reaction_set_basisVector_lb,temp_lb'];
              reaction_set_basisVector_ub=[reaction_set_basisVector_ub,temp_ub'];
              
              single_rev=find(single_reaction<=length(temp_rev));
              if ~isempty(single_rev)
                single_rev=single_reaction(single_rev);
                rev_reaction=Ecoli_Model_specific.S(:,reaction_set_index(single_rev))*-1;
                single_reaction=setdiff(single_reaction,single_rev);
                temp_reactionset=[temp_reactionset;rev_reaction'];
                for k=1:length(single_rev)
                  reaction_basisVector_index{index}=reaction_set_index(single_rev(k));
                  reaction_basisVector_coff{index}=temp_basisVector(single_rev(k))*-1/abs(temp_basisVector(single_rev(k)));
                  index=index+1;
                end
              end
              for_reaction=Ecoli_Model_specific.S(:,reaction_set_index(single_reaction));
              temp_reactionset=[temp_reactionset;for_reaction'];
              for k=1:length(single_reaction)
                  reaction_basisVector_index{index}=reaction_set_index(single_reaction(k));
                  reaction_basisVector_coff{index}=temp_basisVector(single_reaction(k))/abs(temp_basisVector(single_reaction(k)));
                  index=index+1;
              end
              continue;
          end
          temp_reaction_index=find(temp_basisVector ~= 0);
          temp_reaction=[];
          if length(temp_reaction_index)==1
              for k=1:length(temp_reaction_index)
                  if repeat_hash.isEmpty || ~repeat_hash.containsKey(reaction_set_index(temp_reaction_index(k)))
                      if temp_reaction_index(k)<=length(temp_rev)
                          repeat_hash.put(reaction_set_index(temp_reaction_index(k)),-1);
                          temp_reaction=[temp_reaction,temp_reaction_index(k)];
                      else
                          repeat_hash.put(reaction_set_index(temp_reaction_index(k)),1);
                          temp_reaction=[temp_reaction,temp_reaction_index(k)];
                      end
                  else
                      if temp_reaction_index(k)<=length(temp_rev) && repeat_hash.get(reaction_set_index(temp_reaction_index(k)))==1
                          temp_reaction=[temp_reaction,temp_reaction_index(k)];
                          repeat_hash.put(reaction_set_index(temp_reaction_index(k)),0);
                      elseif temp_reaction_index(k)> length(temp_rev) && repeat_hash.get(reaction_set_index(temp_reaction_index(k)))==-1
                          temp_reaction=[temp_reaction,temp_reaction_index(k)];
                          repeat_hash.put(reaction_set_index(temp_reaction_index(k)),0);
                      end
                  end
              end
              if isempty(temp_reaction)
                  continue;
              end
          end
          original_temp_reaction_index=temp_reaction_index;
          temp_reaction_index=reaction_set_index(temp_reaction_index);
          reaction_basisVector_index{index}=temp_reaction_index;
          
          single_rev=find(original_temp_reaction_index<=length(temp_rev));
          if ~isempty(single_rev)
              single_rev=original_temp_reaction_index(single_rev);
              temp_basisVector(single_rev)=temp_basisVector(single_rev)*-1;
          end
          reaction_basisVector_coff{index}=temp_basisVector(original_temp_reaction_index);
          index=index+1;
          
          temp_reaction=zeros(size(Ecoli_Model_specific.S,1),1);
          temp_reaction(metabolite_index)=temp_metabolites;
          temp_reactionset=[temp_reactionset;temp_reaction'];
          temp_index=find(temp_basisVector ~= 0);
          temp_lb=max(sub_model(i).lb(temp_index)./abs(temp_basisVector(temp_index)));
          temp_ub=min(sub_model(i).ub(temp_index)./abs(temp_basisVector(temp_index)));
          reaction_set_basisVector_lb=[reaction_set_basisVector_lb,temp_lb];
          reaction_set_basisVector_ub=[reaction_set_basisVector_ub,temp_ub];
      end
      reaction_set_basisVector{i}=temp_reactionset';
  end
  
  linear_matrix=Ecoli_Model_specific.S(:,noClusterReactionSet);
  linear_matrix=linear_matrix';
  for i=1:length(reaction_set_basisVector)
      linear_matrix=[linear_matrix;1*reaction_set_basisVector{i}'];
  end
  linear_matrix=linear_matrix';
  
  zero_metabolite=[];
  for i=1:size(linear_matrix,1)
      if sum(abs(linear_matrix(i,:)))==0
          zero_metabolite=[zero_metabolite,i];
      end
  end
  metabolite_index=1:size(linear_matrix,1);
  nonzero_metabolite_index=setdiff(metabolite_index,zero_metabolite);
  linear_matrix=linear_matrix(nonzero_metabolite_index,:);
  zero_metaboliteSet={};
  for i=1:length(zero_metabolite)
      zero_metaboliteSet{i}=Ecoli_Model_specific.mets{zero_metabolite(i)};
  end
  
  % reconstruct new model
  for i=1:length(noClusterReactionSet)
      linear_model.rxns{i}=Ecoli_Model_specific.rxns{noClusterReactionSet(i)};
      linear_model.rxns{i}=Ecoli_Model_specific.rxnNames{noClusterReactionSet(i)};
  end
  for i=length(noClusterReactionSet)+1:size(linear_matrix,2)
      linear_model.rxns{i}=strcat('basis_reaction_',num2str(i));
      linear_model.rxnNames{i}=strcat('basis_reaction_',num2str(i));
  end
  for i=1:length(nonzero_metabolite_index)
      linear_model.mets{i}=Ecoli_Model_specific.mets{nonzero_metabolite_index(i)};
      linear_model.metNames{i}=Ecoli_Model_specific.metNames{nonzero_metabolite_index(i)};
  end
  linear_model.S=linear_matrix;
  linear_model.basisVectorIndex=reaction_basisVector_index;
  linear_model.reaction_basisVector_coff=reaction_basisVector_coff;
  linear_model.rev=Ecoli_Model_specific.rev(noClusterReactionSet);
  linear_model.rev(length(noClusterReactionSet)+1:size(linear_matrix,2))=0;
  linear_model.lb=Ecoli_Model_specific.lb(noClusterReactionSet);
  linear_model.lb(length(noClusterReactionSet)+1:size(linear_matrix,2))=reaction_set_basisVector_lb;
  linear_model.ub=Ecoli_Model_specific.ub(noClusterReactionSet);
  linear_model.ub(length(noClusterReactionSet)+1:size(linear_matrix,2))=reaction_set_basisVector_ub;
  linear_model.c=Ecoli_Model_specific.c(noClusterReactionSet);
  linear_model.c(length(noClusterReactionSet)+1:size(linear_matrix,2))=0;
  biomass_index=0;
  for i=1:length(linear_model.rxns)
      if ~isempty(regexp(linear_model.rxns{i},'biomass\sobjective\sfunction','match'))
          biomass_index=i;
          break;
      end
  end
  linear_model.c(biomass_index)=1;
  linear_model.minPathwayStart=length(noClusterReactionSet);
  linear_model.reserveReac=noClusterReactionSet;
  linear_model.minPathwayNum=length(linear_model.basisVectorIndex);
%   linear_model.c(1585)=1;
  linear_model.description='linear model';
  linear_model.b=Ecoli_Model_specific.b(nonzero_metabolite_index);
%   for i=1:length(Nsnp)
%       linear_model.basisvector_set{i}=Nsnp{i};
%       
%       linear_model.basisvector_reaction_lb{i}=sub_model(i).lb;
%       linear_model.basisvector_reaction_ub{i}=sub_model(i).ub;
%   end
new_network=ConnectGraphConstruct(linear_model.S,linear_model.rxns,linear_model.mets,cofactor_hash,linear_model.rev);
fid=fopen('new_netwrok.txt','w');
for i=1:length(new_network.metabolite_reaction_connect)
    temp=new_network.metabolite_reaction_connect{i};
    if ~isempty(temp)
        for j=1:length(temp)
            fprintf(fid,'%d\t%d\n',i,temp(j));
        end
    else
%         fprintf(fid,'%d\t%d\n',i,i);
    end
end
fclose(fid);
%% linear network connection supply
  
% get the infeasible metabolites set
old_linear_model=linear_model;
update_model=linear_model;
test_meta=find(update_model.S(:,biomass_index)<0);
% test_meta=find(update_model.S(:,517)~=0);
orginal_coffi=update_model.S(test_meta,biomass_index);
infeasible_input=[];
undate_len=0;
while(1)
     for i=1:length(test_meta)
          solution=optimizeCbModel(update_model,[],[],1);
          
          if solution.f~=0
              if i>1
                disp(solution.f)
                disp(solution.x(1391))
                infeasible_input=[infeasible_input,test_meta(i-1)];
              end
              break;
          else
              update_model.S(test_meta(i),biomass_index)=0;
          end
     end
     if length(infeasible_input)==undate_len
         break;
     end
     update_model.S(test_meta,biomass_index)=orginal_coffi;
     update_model.S(infeasible_input,biomass_index)=0;
     undate_len=length(infeasible_input);
end


% infeasible set regulization
update_model=linear_model;
coeffi=update_model.S(infeasible_input,biomass_index);
candidate_index=1:size(update_model.S,1);
candidate_index=setdiff(candidate_index,infeasible_input);
validate_index=zeros(length(candidate_index),length(infeasible_input));
flag_all=0;
infeasible_process=infeasible_input;

for i=1:length(candidate_index)
    update_model=linear_model;
    temp=zeros(1,size(update_model.S,1));
    temp(candidate_index(i))=1;
    update_model.S=[update_model.S';temp];
    update_model.S=update_model.S';
    update_model.ub=[update_model.ub;1000];
    update_model.lb=[update_model.lb;0];
    update_model.c=[update_model.c;0];
    update_model.c(biomass_index)=1;
    update_model.rev=[update_model.rev;0];
    
    for j=1:length(infeasible_input)
        update_model.S(infeasible_input,biomass_index)=0;
        update_model.S(infeasible_input(j),biomass_index)=coeffi(j);
        test_solution=optimizeCbModel(update_model,[],[],1);
        validate_index(i,j)=test_solution.f;
    end
end


new_meta_set=[];
for i=1:length(infeasible_input)
    temp=max(validate_index(:,i));
%     temp=0.95*temp;
    if temp<1
        new_meta_set=[new_meta_set,infeasible_input(i)];
    else
        sub_temp=find(validate_index(:,i)>=temp);
        sub_temp=candidate_index(sub_temp);
        new_meta_set=[new_meta_set,sub_temp(1)];
%         new_meta_set=[new_meta_set,sub_temp];
    end
end


new_meta_set=unique(new_meta_set);
validate_index=zeros(length(new_meta_set),size(linear_model.S,1));
for i=1:length(new_meta_set)
    validate_index(i,new_meta_set(i))=1;
end

linear_model.S=[linear_model.S';validate_index];
linear_model.S=linear_model.S';
linear_model.lb=[linear_model.lb;zeros(size(validate_index,1),1)];
linear_model.ub=[linear_model.ub;1000*ones(size(validate_index,1),1)];
linear_model.c=[linear_model.c;zeros(size(validate_index,1),1)];
linear_model.rev=[linear_model.rev;zeros(size(validate_index,1),1)];
for i=length(linear_model.rxns)+1:size(linear_model.S,2)
    linear_model.rxns{i}=strcat('basis_reaction_',num2str(i));
    linear_model.rxnNames{i}=strcat('basis_reaction_',num2str(i));
end
  %% the flux balance analysis for the linear model by the MILP solver
  solution=optimizeCbModel(linear_model,[],[],1);
  Ecoli_Model.c(1005)=1;
  solution1=optimizeCbModel(Ecoli_Model);
  
  linear_FBA_acReac_num=find(solution.x~=0);
  linear_FBA_acReac_index=[];
  linear_FBA_acReac_flux=[];
  for i=1:length(linear_FBA_acReac_num)
      if linear_FBA_acReac_num(i)<=linear_model.minPathwayStart
          linear_FBA_acReac_index=[linear_FBA_acReac_index,linear_model.reserveReac(linear_FBA_acReac_num(i))];
          linear_FBA_acReac_flux=[linear_FBA_acReac_flux,solution.x(linear_FBA_acReac_num(i))];
      elseif linear_FBA_acReac_num(i)-linear_model.minPathwayStart <= linear_model.minPathwayNum
          linear_FBA_acReac_index=[linear_FBA_acReac_index,linear_model.basisVectorIndex{linear_FBA_acReac_num(i)-linear_model.minPathwayStart}];
          linear_FBA_acReac_flux=[linear_FBA_acReac_flux,[solution.x(linear_FBA_acReac_num(i))*linear_model.reaction_basisVector_coff{linear_FBA_acReac_num(i)-linear_model.minPathwayStart}]'];
      else
          linear_FBA_acReac_index=[linear_FBA_acReac_index,linear_FBA_acReac_num(i)];
          linear_FBA_acReac_flux=[linear_FBA_acReac_flux,solution.x(linear_FBA_acReac_num(i))];
      end
  end
  FBA_acReac_index=find(solution1.x~=0);
  unique_linear_FBA=unique(linear_FBA_acReac_index);
  unique_linear_flux=[];
  for i=1:length(unique_linear_FBA)
      index=find(linear_FBA_acReac_index==unique_linear_FBA(i));
      unique_linear_flux=[unique_linear_flux,sum(linear_FBA_acReac_flux(index))];
  end
  unique_linear_index=find(unique_linear_FBA<=length(specific_reaction));
  unique_linear_FBA=unique_linear_FBA(unique_linear_index);
  unique_linear_flux=unique_linear_flux(unique_linear_index);
  unique_linear_FBA=specific_reaction(unique_linear_FBA);
  
  recontribut_model.S=Ecoli_Model.S;
  recontribut_model.c=Ecoli_Model.c;
  recontribut_model.b=Ecoli_Model.b;
  recontribut_model.lb=Ecoli_Model.lb;
  recontribut_model.lb(unique_linear_FBA)=unique_linear_flux-100;
%   recontribut_model.lb(FBA_acReac_index)=solution1.x(FBA_acReac_index)-0.5;
  recontribut_model.lb(1005)=0;
  recontribut_model.ub(1235)=0;
  recontribut_model.ub(1233)=0;
%     recontribut_model.lb(2377)=0;
    recontribut_model.lb(find(recontribut_model.lb<-1000))=-1000;
  recontribut_model.ub=Ecoli_Model.ub;
  recontribut_model.ub(unique_linear_FBA)=unique_linear_flux+100;
%   recontribut_model.ub(FBA_acReac_index)=solution1.x(FBA_acReac_index)+0.5;
  recontribut_model.ub(1005)=1000;
%   recontribut_model.ub(2377)=1000;
%   recontribut_model.ub(2369)=1000;
%   recontribut_model.ub(2372)=1000;
  recontribut_model.ub(1235)=1000;
  recontribut_model.ub(1233)=1000;
  recontribut_model.ub(find(recontribut_model.ub>1000))=1000;
  solutionr=optimizeCbModel(recontribut_model,[],[],1);
 %% C13 flux analysis
 reaction_hash=java.util.Hashtable;
 for i=1:length(Ecoli_Model.rxns)
     reaction_hash.put(Ecoli_Model.rxns{i},i);
 end
 file_id=fopen('D:\fluxmode\flux analysis\C13\fomulate_new.txt');
 index=1;
 C13_data_index=[];
 while(~feof(file_id))
     line=fgetl(file_id);
     str=regexp(line,'\t','split');
     if strcmp(str{1},'reactionname') || length(str)<5
         continue;
     end
     sub_str=str{1};
     sub_str=sub_str(3:end);
     disp(sub_str)
     if length(sub_str)<2
         continue;
     end
     if reaction_hash.containsKey(sub_str)
         C13_data_index(index,1)=reaction_hash.get(sub_str);
         temp=[];
         for i=2:length(str)
             if ~isempty(str{i})
                 temp=[temp,str2num(str{i})];
             else
                 temp=[temp,0];
             end
         end
         C13_data_index(index,2:5)=temp;
         index=index+1;
     end
 end
 fclose(file_id);
 
 C13data_OnOff=C13_data_index(:,2);
 C13data_OnOff_label=C13data_OnOff;
 C13data_OnOff_label(find(C13data_OnOff_label~=0))=1;
 FBAdata_OnOff=solution1.x(C13_data_index(:,1));
 FBAdata_OnOff_label=FBAdata_OnOff;
 FBAdata_OnOff_label(find(FBAdata_OnOff_label~=0))=1;
 sum_label=FBAdata_OnOff_label+C13data_OnOff_label;
 length(find(sum_label==2))
 corr(abs(FBAdata_OnOff(find(sum_label==2))),abs(C13data_OnOff(find(sum_label==2))))
 corr(abs(FBAdata_OnOff(find(sum_label>0))),abs(C13data_OnOff(find(sum_label>0))))
 corr(log(abs(FBAdata_OnOff(find(sum_label>0)))+1),log(abs(C13data_OnOff(find(sum_label>0)))+1))
 
 figure(1)
 para1=plot(log(abs(FBAdata_OnOff(find(sum_label>0)))+1),log(abs(C13data_OnOff(find(sum_label>0)))+1),'b+',0:0.5:7,0:0.5:7,'--k');
 axis([0 7 0 7]);
 set(gca,'Fontsize',15,'Fontname','Timesnewroman','FontWeight','bold','linew',2);
 set(para1(1),'LineWidth',2);
 set(para1(2),'LineWidth',2);
 MSE_FBA=sqrt((abs(FBAdata_OnOff(find(sum_label>0)))-abs(C13data_OnOff(find(sum_label>0))))'*(abs(FBAdata_OnOff(find(sum_label>0)))-abs(C13data_OnOff(find(sum_label>0)))));
 MSE_FBA=MSE_FBA/length(find(sum_label>0))
 
 C13data_OnOff=C13_data_index(:,2);
 C13data_OnOff_label=C13data_OnOff;
 C13data_OnOff_label(find(C13data_OnOff_label~=0))=1;
 FBAdata_OnOff=solutionr.x(C13_data_index(:,1));
 FBAdata_OnOff_label=FBAdata_OnOff;
 FBAdata_OnOff_label(find(FBAdata_OnOff_label~=0))=1;
 sum_label=FBAdata_OnOff_label+C13data_OnOff_label;
 length(find(sum_label==2))
 corr(abs(FBAdata_OnOff(find(sum_label==2))),abs(C13data_OnOff(find(sum_label==2))))
 corr(abs(FBAdata_OnOff(find(sum_label>0))),abs(C13data_OnOff(find(sum_label>0))))
 corr(log(abs(FBAdata_OnOff(find(sum_label>0)))+1),log(abs(C13data_OnOff(find(sum_label>0)))+1))
 figure(2)
 para1=plot(log(abs(FBAdata_OnOff(find(sum_label>0)))+1),log(abs(C13data_OnOff(find(sum_label>0)))+1),'r+',0:0.5:7,0:0.5:7,'--k');
 axis([0 7 0 7]);
 set(gca,'Fontsize',15,'Fontname','Timesnewroman','FontWeight','bold','linew',2);
 set(para1(1),'LineWidth',2);
 set(para1(2),'LineWidth',2);
 MSE_LFBA=sqrt((abs(FBAdata_OnOff(find(sum_label>0)))-abs(C13data_OnOff(find(sum_label>0))))'*(abs(FBAdata_OnOff(find(sum_label>0)))-abs(C13data_OnOff(find(sum_label>0)))));
 MSE_LFBA=MSE_LFBA/length(find(sum_label>0))
 %% end C13
 %% C13 flux confidence range
 Flux_model=Ecoli_Model;
 bound_LFBA=[recontribut_model.lb';recontribut_model.ub']';
 solution_FBA=solution1;
 solution_LFBA=solutionr;
end