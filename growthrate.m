gene_growthrat_hash = java.util.Hashtable;
input_nutrient_file=fopen('.\data\growth_rate.txt');
while(~feof(input_nutrient_file))
    line=fgetl(input_nutrient_file);
    if(~isempty(line))
        str=regexp(line,'\t','split');
        if length(str)<2
            continue;
        end
        gene_growthrat_hash.put(str{1},str2num(str{2}));
    end
end
fclose(input_nutrient_file);

[Flux_model,bound_LFBA,solution_FBA,solution_LFBA]=linear_FBA_function();
Flux_LFBA_model=Flux_model;
Flux_LFBA_model.lb=bound_LFBA(:,1);
Flux_LFBA_model.ub=bound_LFBA(:,2);
FBA_KO_growthrate=[];
MOMA_KO_growthrate=[];
LFBA_KO_growthrate=[];
LFBA_MOMA_KO_growthrate=[];
measure_growthrate=[];
genelist={};
index=1;
for i=1:length(Flux_model.genes)
    if gene_growthrat_hash.containsKey(Flux_model.genes{i})
        growthRate=gene_growthrat_hash.get(Flux_model.genes{i});
        measure_growthrate=[measure_growthrate,growthRate];
        genelist{index}=Flux_model.genes{i};
        index=index+1;
    end
end
 [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(Flux_model, 'FBA', genelist, solution_FBA);
 FBA_KO_growthrate=[grRatio';grRateKO';ones(1,length(grRateKO))*grRateWT']';
 [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(Flux_model, 'lMOMA', genelist, solution_FBA);
 MOMA_KO_growthrate=[grRatio';grRateKO';ones(1,length(grRateKO))*grRateWT']';
%  [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(Flux_LFBA_model, 'FBA', genelist);
 [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(Flux_LFBA_model, 'FBA', genelist, solution_LFBA);
 LFBA_KO_growthrate=[grRatio';grRateKO';ones(1,length(grRateKO))*grRateWT']';
 [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(Flux_LFBA_model, 'lMOMA', genelist, solution_LFBA);
 LFBA_MOMA_KO_growthrate=[grRatio';grRateKO';ones(1,length(grRateKO))*grRateWT']';
corr(measure_growthrate',FBA_KO_growthrate(:,1))
corr(measure_growthrate',MOMA_KO_growthrate(:,1))
index_temp=find(isnan(LFBA_KO_growthrate(:,1))==0);
data_temp=LFBA_KO_growthrate(:,1);
corr(measure_growthrate(index_temp)',data_temp(index_temp))
% corr(measure_growthrate',LFBA_KO_growthrate(:,1))
index_temp=find(isnan(LFBA_MOMA_KO_growthrate(:,1))==0);
data_temp=LFBA_MOMA_KO_growthrate(:,1);
corr(measure_growthrate(index_temp)',data_temp(index_temp))
% corr(measure_growthrate',LFBA_MOMA_KO_growthrate(:,1))
 figure(3)
 lb_axis=min(min([measure_growthrate',FBA_KO_growthrate(:,1)]));
 ub_axis=max(max([measure_growthrate',FBA_KO_growthrate(:,1)]));
 para1=plot(measure_growthrate',FBA_KO_growthrate(:,1),'r+',lb_axis:0.1:ub_axis,lb_axis:0.1:ub_axis,'--k');
 axis([lb_axis ub_axis lb_axis ub_axis]);
 set(gca,'Fontsize',15,'Fontname','Timesnewroman','FontWeight','bold','linew',2);
 set(para1(1),'LineWidth',2);
 set(para1(2),'LineWidth',2);
 figure(4)
 lb_axis=min(min([measure_growthrate',MOMA_KO_growthrate(:,1)]));
 ub_axis=max(max([measure_growthrate',MOMA_KO_growthrate(:,1)]));
 para1=plot(measure_growthrate',MOMA_KO_growthrate(:,1),'r+',lb_axis:0.1:ub_axis,lb_axis:0.1:ub_axis,'--k');
 axis([lb_axis ub_axis lb_axis ub_axis]);
 set(gca,'Fontsize',15,'Fontname','Timesnewroman','FontWeight','bold','linew',2);
 set(para1(1),'LineWidth',2);
 set(para1(2),'LineWidth',2);
 figure(5)
 lb_axis=min(min([measure_growthrate',LFBA_KO_growthrate(:,1)]));
 ub_axis=max(max([measure_growthrate',LFBA_KO_growthrate(:,1)]));
 para1=plot(measure_growthrate',LFBA_KO_growthrate(:,1),'r+',lb_axis:0.1:ub_axis,lb_axis:0.1:ub_axis,'--k');
 axis([lb_axis ub_axis lb_axis ub_axis]);
 set(gca,'Fontsize',15,'Fontname','Timesnewroman','FontWeight','bold','linew',2);
 set(para1(1),'LineWidth',2);
 set(para1(2),'LineWidth',2);
 figure(6)
 lb_axis=min(min([measure_growthrate',LFBA_MOMA_KO_growthrate(:,1)]));
 ub_axis=max(max([measure_growthrate',LFBA_MOMA_KO_growthrate(:,1)]));
 para1=plot(measure_growthrate',LFBA_MOMA_KO_growthrate(:,1),'r+',lb_axis:0.1:ub_axis,lb_axis:0.1:ub_axis,'--k');
 axis([lb_axis ub_axis lb_axis ub_axis]);
 set(gca,'Fontsize',15,'Fontname','Timesnewroman','FontWeight','bold','linew',2);
 set(para1(1),'LineWidth',2);
 set(para1(2),'LineWidth',2);