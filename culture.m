function [ind_shares,trait_shares,ind_skew,trait_skew,traits_time_converged,inds_time_converged ...
    ,two_inds_time_converged,indicators_interacted_avg,good_interactions_avg, ...
    interactions_avg,unique_indicators,unique_traits,fitness_skew, ...
    swings,avg_ind_entropy,avg_trait_entropy,trait_modularity,ind_modularity,varargout] ...
    =culture(config,range,numagents,numgenerations,extra,r,plots,adjmat,diffs)


%CONFIG is a number from 0 to 512 that translates to binary where each position
%represents presence or absence of  certain features
%1st bit: 0=no 'learning'/1='learning'
%2nd: Naive/ACT memory        
%3rd: no space/space
%4th: singe/multiple indicators
%5th: single/multiple prefs
%6th: single/multiple traits
%7th: b
%8th: no network/ network
%9th: lamarckian  preferences/ genetic preferences
% Most of the bits were fixed throghout experimentation for Revay & Cioffi
% (2018). The following configurations were tested: 010000110, 010000010,
% 010000011.
% Furthermore, many of the bits are deprecated, i.e. they do not affect the
%anymore in any way, and are simply remnants from older versions of the
%model
%OTHER PARAMETERS:
%RANGE: number of possible traits and tags 
%NUMAGENTS: number of agents in the simulation
%NUMGENERATIONS: number of generations in the simulation
%EXTRA: in the ACT-R case this refers to the BLA threshold. In the naive
%case this refers to the memory 'length'.
%R: ratio of interaction/evolution (horizontal to vertical transmission)
%PLOTS: if 'plots' then charts are plotted at the end of the simulations
%ADJMAT: adjacency matrix pf the network on which the population lives
%DIFFS: a matrix derived from the ADJMAT which gives the difference between
%the expected and actual link counts. Used for modularity calculations.
%%% The model outputs a large number of variables, however, there are
%%% other variables on top of thes that the model keeps track of, that the 
%%%user can output if they wish.






%initialize these as the sizes for matrices holding info on the
%indicators/traits/prefernces
if strcmp(config(4),'1')
 ni=2;
else
 ni=1;  
end

if strcmp(config(6),'1')
 nt=2;
else
 nt=1;
end

if strcmp(config(5),'1')
 np=2;
else
 np=1;
end

%these are all book-keeping vector for in-run and output stats:
trait_modularity=nan(numgenerations,1);
ind_modularity=nan(numgenerations,1);
two_m=length(find(adjmat));
ind_changes=[];
trait_changes=[];
neighborhood=[];
indicator_entropy=nan(range,numgenerations,ni,nt);
trait_entropy=nan(range,numgenerations,nt,ni);
local_entropy=zeros(numgenerations,numagents,ni,nt);
indicators_interacted_avg=zeros(numgenerations,ni);
good_interactions_avg=zeros(numgenerations,1);
interactions_avg=zeros(numgenerations,1);
unique_indicators=zeros(numgenerations,ni);
unique_traits=zeros(numgenerations,nt);
fitness_skew=zeros(numgenerations,1);
swings=0;
ind_counts=zeros(numgenerations,range,ni);
trait_counts=zeros(numgenerations,range,nt);
pref_counts=zeros(numgenerations,range,np);
ind_shares=zeros(ni,numgenerations);
trait_shares=zeros(nt,numgenerations);
ind_skew=zeros(ni,numgenerations);
trait_skew=zeros(nt,numgenerations);
pref_shares=zeros(np,numgenerations);
traits_converged=0;
traits_time_converged=2001;
inds_converged=0;
inds_time_converged=2001;
two_ind_shares=zeros(1,numgenerations);
two_inds_time_converged=2001;
two_inds_converged=0;
best_trait_old=0;
best_trait_new=0;
shares_per_inds=zeros(numgenerations,range,ni,nt);
trait_counts_per_inds=cell(range,ni);
for i=1:range
    for k=1:ni
        trait_counts_per_inds{i,k}=zeros(numgenerations,range,nt);
    end
end


%%initialize agent variables:
indicators=randi(range,numagents,ni);
traits=randi(range,numagents,nt);
fitness=zeros(numagents,1);
good_interactions=zeros(numagents,1);
interactions=zeros(numagents,1);
indicators_interacted=zeros(numagents,range,ni);

%initialize the agent population as a structure from the previously
%initiliazed variables:
%minimal agent genotype:
agents=struct('fitness',fitness,'interactions',interactions,...
    'good_interactions',good_interactions,'indicators',indicators,...
    'traits',traits,'indicators_interacted',indicators_interacted);

%%% initialize the old 'learning' variables for th agents
if strcmp(config(1),'1') 
    learnfreq_avg=zeros(numgenerations,1);
    learnfreq=rand(numagents,1);
    agents.learnfreq=learnfreq;
    if strcmp(config(2),'0')
        learned_traits=randi(range,numagents,nt,extra);
        agents.learned_traits=learned_traits;  
    else
        learned_n=zeros(numagents,range,nt);
        learned_t=zeros(numagents,range,nt);
        agents.learned_n=learned_n;
        agents.learned_t=learned_t;
    end
end

%initialize the variables needed for ACT-R memory:
if strcmp(config(2),'1')
    pref_n=randi(20,numagents,range,np); %counts of good interactions
    pref_t=randi(20*numagents,numagents,range,np); %times of good interactions
    bad_n=randi(20,numagents,range,np); %counts of bad interactions
    bad_t=randi(20*numagents,numagents,range,np);%times of bad interactions 
    pref_bla_avgs=zeros(100,range,np);
     agents.pref_n=pref_n;
     agents.pref_t=pref_t;
     agents.bad_n=bad_n;
     agents.bad_t=bad_t;
%otherwise initialize variables for 'naive' memory
else
    prefs=randi(range,numagents,extra,np); %tags from past good memories
    preflens=randi(floor(r*numagents),numagents,extra,np); %times they have been in memory
    agents.prefs=prefs;
    agents.preflens=preflens;
end

%initiliaze neighoborrhoods and book-keeping variables for the moore
%neighborhood version of the model:
if strcmp(config(3),'1')
   neighborhood=make_moore_nbhd(numagents);
    num_ind_clusters=zeros(numgenerations/10,ni);
    num_trait_clusters=zeros(numgenerations/10,nt);
    ind_cluster_sizes=cell(numgenerations/10,ni);
    trait_cluster_sizes=cell(numgenerations/10,nt);
    trait_changes=zeros(numagents,nt);
    ind_changes=zeros(numagents,ni);
    agents.nbhd=neighborhood;
  
end

% extra book-keeping variables for the multiple-tag scenario
if strcmp(config(4),'1')
    inds_time_converged=[2001 2001];
    inds_converged=[0 0];
    two_ind_shares=zeros(2,numgenerations);
    two_inds_time_converged=[2001 2001];
    two_inds_converged=[0 0];
    if strcmp(config(5),'1')
       p2i_map=randi(2,numagents,2); 
       agents.p2i=p2i_map;
    end
    
end

%extra book-keeping variables for the multiple preferences scenario 
%(THIS MIGHT BE DEPRECATED)
% if strcmp(config(5),'1')
%    if strcmp(config(2),'1')
%      pref_n=zeros(numagents,range,2);
%      pref_t=inf(numagents,range,2);
%      bad_n=zeros(numagents,range,2);
%      bad_t=inf(numagents,range,2);
%      pref_bla_avgs=zeros(numgenerations,range,2);
%      agents.pref_n=pref_n;
%      agents.pref_t=pref_t;
%      agents.bad_n=bad_n;
%      agents.bad_t=bad_t;
%    else
%      prefs=randi(range,numagents,extra,2);
%      agents.prefs=prefs;
%    end
% end

%extra book-keeping for the multiple trait case:
if strcmp(config(6),'1')
    traits_time_converged=[2001 2001];
    traits_converged=[0 0];
    best_trait_old=[0 0];
    best_trait_new=[0 0];
    swings=[0 0];
end

%visualization details in case the user wants to plot results:
if strcmp(plots,'plots')
    % jsut some neat colormaps for the visualizatoins:
                trait_map=[ 0    0.0794    0.9603;
         0    0.1746    0.9127;
         0    0.2698    0.8651;
         0    0.3651    0.8175;
         0    0.4603    0.7698;
         0    0.5556    0.7222;
         0    0.6508    0.6746;
         0    0.7460    0.6270;
         0    0.8413    0.5794;
         0    0.9365    0.5317];
     
    ind_map=[1.0000    0.0794    0.9206;
    1.0000    0.1746    0.8254;
    1.0000    0.2698    0.7302;
    1.0000    0.3651    0.6349;
    1.0000    0.4603    0.5397;
    1.0000    0.5556    0.4444;
    1.0000    0.6508    0.3492;
    1.0000    0.7460    0.2540;
    1.0000    0.8413    0.1587;
    1.0000    0.9365    0.0635];

    pref_map=[
      1 1 0;
      0.84 0.84 0;
      0.68 0.68 0;
      1 0 1;
      0.84 0 0.84;
      0.68 0 0.68;
      0 1 1;
      0 0.84 0.84;
      0 0.68 0.68;
      1 0 0;
      0.84 0 0;
      0.68 0 0;
      0 1 0;
      0 0.84 0;
      0 0.68 0;
      0 0 1;
      0 0 0.84;
      0 0 0.68;
      0.52 0.52 0;
      0.36 0.36 0;
      0.2 0.2 0;
      0.52 0 0.52;
      0.36 0 0.36;
      0.2 0 0.2;
      0 0.52 0.52;
      0 0.36 0.36;  
      0 0.2 0.2;
      0.36 0.36 0.36; 
      0.2 0.2 0.2;
      0 0 0;
    ];
    agents.prefmapping=nan(numagents,1);
end


%bits 7,8, and 9 do not require any special preallocation




%times at which the network data is written out:
gens=[5 10 20 30 40 50];
%begin generation:
for n=1:numgenerations

   sprintf('generation: %d', n)
   
  % write out network data at the chosen times: 
  if ~isempty(find(gens==n,1))
   fid = fopen(sprintf('sw-gen%d-nodes.csv',n),'w');
   fprintf(fid,'Id,Trait,Indicator\n');
   fclose(fid);
   ids=1:1:numagents;
   nodes=ids';
   nodes=[nodes agents.traits agents.indicators];
   dlmwrite(sprintf('sw-gen%d-nodes.csv',n),nodes,'delimiter',',','-append'); %nodelist (agents with their traits and tags)
   list=edgelist(adjmat,numagents,'unweighted','zeros'); 
   csvwrite(sprintf('sw-edgelist-gen%d',n),list); %edgelist (links between agents)
  end
   
   % visualization in the grid scenario:
   if strcmp(config(3),'1')
       
       trait_grid=reshape(agents.traits(:,1),sqrt(numagents),sqrt(numagents));
       if strcmp(plots,'plots')
        
        %draw the trait grid:
        figure(1)
        image(trait_grid);
        colormap(trait_map)
        colorbar
        drawnow;
               
        %draw the tag grid: 
        ind_grid=reshape(agents.indicators(:,1),sqrt(numagents),sqrt(numagents));
        figure(2)
        image(ind_grid);
        colormap(ind_map)
        colorbar
        drawnow;
        
        %draw the preference grid:
        pref_grid=reshape(agents.prefmapping(:,1),sqrt(numagents),sqrt(numagents));
        figure(3)
        image(pref_grid);
        colormap(pref_map)
        colorbar
        drawnow; 
       end
       
   end
       

   %each round:
   for j=1:numagents*r
      
       i=randi(numagents); %select random agents
       
       if strcmp(config(6),'1')
          trait_locus=randi(2,1); %randomly select cultural trait
       else
          trait_locus=1;
       end
       
       if strcmp(config(5),'1') && strcmp(config(6),'1')
          pref_locus=trait_locus; %if two sets of preferences and traits, ...
          %get the one associated with the selected trait locus       
       elseif strcmp(config(5),'1') && strcmp(config(6),'0')
          pref_locus=randi(2,1); %if two sets of preferneces but only one ...
          %trait, choose the preference locus randomly
       else 
          pref_locus=1;
       end    
       
       
       if strcmp(config(5),'1') && strcmp(config(4),'1')
          ind_locus=agents.p2i(i,pref_locus); %if two sets of tags and ...
          %preferences assign the tag locus to the preference locuse based...
          %on the agent's mapping
       else    
          ind_locus=1;
       end
       
       %old learning mechanism, not used in experimentation: 
       %agents=learn(i,agents,config,pref_locus,ind_locus,extra,extra2,n,range,nt);
       
       %have the agent select potential partners based on their
       %preferences:
       [potential_partners,preferred_ind_variant,divider]=get_potential_partners(i,agents,config,pref_locus,ind_locus,extra,n,range,adjmat);
       
       %count this as an interaction, if there are any potential partners:
       if ~isempty(potential_partners)
         agents.interactions(i)=agents.interactions(i)+1;
       end
       
       %select an actual partner from the set of candidates:
       [partner,ind_locus,partner_ind_locus,agents]=get_actual_partner(i,agents,potential_partners,ind_locus,pref_locus,config,n,divider);
   
       %if none is selected give the agents -1 to fitness:    
       if isempty(partner)
          agents.fitness(i)=agents.fitness(i)-1;
          continue
       end

       %otherwise interact with the selected agent:
       agents=interact(i,agents,partner,ind_locus,partner_ind_locus,trait_locus,pref_locus,preferred_ind_variant,config);
       
       %and record interaction stats:
       agents.indicators_interacted(i,agents.indicators(partner),ind_locus)=1;
       agents.indicators_interacted(partner,agents.indicators(i),partner_ind_locus)=1;
       if strcmp(config(2),'1') %ACT-R memory
           agents.pref_t=agents.pref_t+1; %tick the time since last activation
           agents.bad_t=agents.bad_t+1;
       end
       
       if(strcmp(config(2),'0')) %naive memory
           agents.preflens=agents.preflens+1; %tick the age of the memorized info
           agents.prefs(agents.preflens>r*numagents)=nan; %scrap info that is too old
       end
   end
  
    %once all the rounds are done, create new generation of agents:
    [offspring,trait_changes,ind_changes]=create_offspring(range,numagents,ni,nt,np,extra,r,agents,ind_changes,trait_changes,neighborhood,config,n,adjmat); 

    if strcmp(plots,'plots') && strcmp(config(3),'1')
        offspring.prefmapping=agents.prefmapping;
    end

    %mark generation level statistics:
    good_interactions_avg(n)=mean(agents.good_interactions);
    interactions_avg(n)=mean(agents.interactions);
    this_indicators_interacted=zeros(numagents,ni);
    for j=1:ni
        for i=1:numagents
           this_indicators_interacted(i,j)=length(find(agents.indicators_interacted(i,:,j)));
        end
        indicators_interacted_avg(n,j)=mean(this_indicators_interacted(:,j));
        unique_indicators(n,j)=length(unique(agents.indicators(:,j)));
    end
    for j=1:nt
        unique_traits(n,j)=length(unique(agents.traits(:,j)));
    end
    
    for j=1:ni
        for i=1:range
            indices=find(agents.indicators(:,j)==i);
            for k=1:nt
                trait_counts_per_inds{i,j}(n,:,k)=histc(agents.traits(indices,k),1:1:range);
                shares_per_inds(n,i,j,k)=max(trait_counts_per_inds{i,j}(n,:,k))/sum(trait_counts_per_inds{i,j}(n,:,k));
            end
        end   
    end
    
    fitness_skew(n)=skewness(agents.fitness);
    
    agents=offspring;
   
    %calculate local entropy of tags in the grid scenario:
    if strcmp(config(3),'1') 
        for i=1:numagents
            for j=1:ni
                for k=1:nt
                    local_inds=agents.indicators(agents.nbhd(i,:),j);
                    local_traits=agents.traits(agents.nbhd(i,:),k);
                    u_inds=unique(local_inds);
                    for m=1:length(u_inds)
                      local_entropy(n,i,j,k)=local_entropy(n,i,j,k)+length(local_inds(local_inds==u_inds(m)))*shannon_entropy(local_traits,local_inds,u_inds(m),range,0);
                    end
                    local_entropy(n,i,j,k)=local_entropy(n,i,j,k)/8;
                end
            end
        end
    end
    
    %calculate local entropy of tags in the network scenario:
     if strcmp(config(8),'1') 
        for i=1:numagents
            for j=1:ni
                
                for k=1:nt
                    
                    local_inds=agents.indicators(adjmat(i,:)==1,j);
                    local_traits=agents.traits(adjmat(i,:)==1,k);
                    u_inds=unique(local_inds);
                    for m=1:length(u_inds)
                      local_entropy(n,i,j,k)=local_entropy(n,i,j,k)+length(local_inds(local_inds==u_inds(m)))*shannon_entropy(local_traits,local_inds,u_inds(m),range,0);
                    end
                    if ~isempty(local_inds)
                        local_entropy(n,i,j,k)=local_entropy(n,i,j,k)/length(local_inds);
                    else
                       local_entropy(n,i,j,k)=0;
                    end
                end
            end
        end
     end   
    
    %keep track of the distributions of tags, traits, and preferences and
    %various derived statistics:
    for j=1:ni
        ind_counts(n,:,j)=histc(agents.indicators(:,j),1:1:range);
        ind_skew(j,n)=skewness(ind_counts(n,:,j));
        ind_counts_now=ind_counts(n,:,j);
        max_ind=max(ind_counts_now);
        ind_counts_now(find(ind_counts_now==max(ind_counts_now),1))=[];
        second_max_ind=max(ind_counts_now);
        ind_shares(j,n)=max_ind/sum(ind_counts(n,:,j));
        two_ind_shares(j,n)=(max_ind+second_max_ind)/sum(ind_counts(n,:,j));
        if ind_shares(j,n)>0.9 && inds_converged(j)==0
            inds_time_converged(j)=n;
            inds_converged(j)=1;
        end 
        if two_ind_shares(j,n)>0.9 && two_inds_converged(j)==0 && ind_shares(j,n)<0.9;
            two_inds_time_converged(j)=n;
            two_inds_converged(j)=1;
            inds_converged(j)=0;
            inds_time_converged(j)=2001;
        elseif two_ind_shares(j,n)<0.9 && two_inds_converged(j)==1;
            two_inds_time_converged(j)=2001;
            two_inds_converged(j)=0;
        elseif ind_shares(j,n)>0.9 && two_inds_converged(j)==1;
            two_inds_time_converged(j)=2001;
            two_inds_converged(j)=0;
        end
    end
    for j=1:nt
        trait_counts(n,:,j)=histc(agents.traits(:,j),1:1:range);
        trait_shares(j,n)=max(trait_counts(n,:,j))/sum(trait_counts(n,:,j));
        if (n==10)
           sprintf('yo') 
        end
        trait_skew(j,n)=skewness(trait_counts(n,:,j));
        best_trait_old(j)=best_trait_new(j);
        [~,best_trait_new(j)]=max(trait_counts(n,:,j));
        if trait_shares(j,n)>0.9 && traits_converged(j)==0
            traits_time_converged(j)=n;
            traits_converged(j)=1;
        elseif traits_converged(j)==1 && best_trait_old(j)~=best_trait_new(j)
            traits_time_converged(j)=2001;
            traits_converged(j)=0;
            swings(j)=swings(j)+1;
        end
    end
    
    for j=1:np
        if strcmp(config(2),'1')
            pref_bla_avgs(n,:,j)=mean(log(2*agents.pref_n(:,j).*agents.pref_t(:,j).^(-0.5)));
            if min(pref_bla_avgs(n,:,j))<0
                pref_bla_avgs_plus=pref_bla_avgs(n,:,j)-min(pref_bla_avgs(n,:,j));
            else 
                pref_bla_avgs_plus=pref_bla_avgs(n,:,j);

            end
            pref_shares(j,n)=max(pref_bla_avgs_plus)/sum(pref_bla_avgs_plus);
        else
            prefs=nan(numagents*extra,1);
            for k=1:extra
               prefs((k-1)*numagents+1:k*numagents)=agents.prefs(:,k,j);
            end


            pref_counts(n,:,j)=histc(prefs,1:1:range);
            pref_shares(j,n)=max(pref_counts(n,:,j))/sum(pref_counts(n,:,j));            
        end
    end
    
   %calculate tag and trait modularity in the grid scenario:
   if strcmp(config(3),'1')

       trait_modularity(n)=clustering2(trait_grid,trait_counts(n,:,1),agents,'grid','traits',diffs);
       ind_modularity(n)=clustering2(ind_grid,ind_counts(n,:,1),agents,'grid','inds',diffs);

   %calculate tag and trait modularity in the network scenario:
   elseif strcmp(config(8),'1')

       coeff=clustering2(adjmat,trait_counts(n,:,1),agents,'graph','traits',diffs);
       trait_modularity(n)=coeff/two_m;
       coeff=clustering2(adjmat,ind_counts(n,:,1),agents,'graph','inds',diffs);
       ind_modularity(n)=coeff/two_m;
   end
    
   %calculate global entropy:
    for i=1:range
        for j=1:ni
            for k=1:nt
                indicator_entropy(i,n,j,k)=shannon_entropy(agents.traits(:,k),agents.indicators(:,j),i,range,10);
                trait_entropy(i,n,k,j)=shannon_entropy(agents.indicators(:,j),agents.traits(:,k),i,range,10);
            end
        end    
    end
    
    inds=permute(ind_counts,[2 1 3]);
    traits=permute(trait_counts,[2 1 3]);
    avg_ind_entropy=zeros(numgenerations,ni,nt);
    avg_trait_entropy=zeros(numgenerations,nt,ni);
    for i=1:numgenerations
        for j=1:ni
            for k=1:nt
                avg_ind_entropy(i,j,k)=sum(inds(~isnan(indicator_entropy(:,i,j,k)),i).*indicator_entropy(~isnan(indicator_entropy(:,i,j,k)),i,j))/numagents;
                avg_trait_entropy(i,k,j)=sum(traits(~isnan(trait_entropy(:,i,k,j)),i).*trait_entropy(~isnan(trait_entropy(:,i,k,j)),i,k))/numagents; 
            end
        end
    end
    
    if strcmp(config(1),'1')
       learnfreq_avg(n)=mean(agents.learnfreq); 
    end
   
    if strcmp(plots,'plots')
        for i=1:numagents
            if strcmp(config(2),'1')
                [~,best_pref]=max(log(2*agents.pref_n(i,:).*agents.pref_t(i,:).^(-0.5)));
                if best_pref<=extra
                    agents.prefmapping(i)=3*(best_pref-1)+1;
                elseif best_pref<=extra+1
                    agents.prefmapping(i)=3*(best_pref-1)+2;
                else
                    agents.prefmapping(i)=3*(best_pref-1)+3;
                end 
            else
                mode_pref=mode(agents.prefs(i,:));
                mode_pref_share=length(find(agents.prefs(i,:)==mode_pref))/extra;
                if mode_pref_share<=0.33
                    agents.prefmapping(i)=3*(mode_pref-1)+1;
                elseif mode_pref_share<=0.67
                    agents.prefmapping(i)=3*(mode_pref-1)+2;
                else
                    agents.prefmapping(i)=3*(mode_pref-1)+3;
                end
            end
        end
    end
   
    
end

%dlmwrite(sprintf('small-world-n%d-r%d-k40-beta5-nodes.csv',numagents,r),nodes,'delimiter',',','-append');

    %calculate various extra measures for the grid scenario:
    if strcmp(config(3),'1')
        avg_local_ind_entropy=squeeze(mean(local_entropy,2));
        trait_stickiness=1-trait_changes./100;
        ind_stickiness=1-ind_changes./100;
        avg_ind_stickiness=mean(ind_stickiness);
        avg_trait_stickiness=mean(trait_stickiness);
        ind_cluster_sizes_skew=zeros(numgenerations/10,ni);
        trait_cluster_sizes_skew=zeros(numgenerations/10,nt);
        for i=1:numgenerations/10
            for j=1:ni
                ind_cluster_sizes_skew(i,j)=skewness(ind_cluster_sizes{i,j});
            end
            for j=1:nt
                trait_cluster_sizes_skew(i,j)=skewness(trait_cluster_sizes{i,j});
            end
        end
        varargout=cell(7,1);
        varargout{1}=avg_local_ind_entropy;
        varargout{2}=avg_trait_stickiness;
        varargout{3}=avg_ind_stickiness;
        varargout{4}=num_ind_clusters;
        varargout{5}=num_trait_clusters;
        varargout{6}=ind_cluster_sizes_skew;
        varargout{7}=trait_cluster_sizes_skew;
        if strcmp(config(1),'1')
             varargout{8}=learnfreq_avg; 
        end
    elseif strcmp(config(1),'1')
        varargout{1}=learnfreq_avg;
    elseif strcmp(config(8),'1')
        avg_local_ind_entropy=squeeze(mean(local_entropy,2));
    end

    
%plot anything the user wishes to plot:
     if strcmp(plots,'plots')
%         figure
%         n=1;
%         for k=1:ni
%             for j=1:nt
%                 subplot(ni*nt,1,n)
%                 n=n+1;
%                 for i=1:range
%                     plot(indicator_entropy(i,:,k,j),'color',ind_map(i,:))
%                     hold on
%                 end
%                 title(sprintf('Inds %d for traits %d',k,j))
%                 xlabel('time')
%             end
%             
%         end
%         legend('1','2','3','4','5','6','7','8','9','10','location','northeast')
% 
%         figure
%         n=1;
%         for k=1:nt
%             for j=1:ni
%                 subplot(ni*nt,1,n)
%                 n=n+1;
%                 for i=1:range
%                     plot(trait_entropy(i,:,k,j),'color',trait_map(i,:))
%                     hold on
%                 end
%                 title(sprintf('Traits %d for inds %d',k,j))
%                 xlabel('time')
%             end
%             
%         end
%         legend('1','2','3','4','5','6','7','8','9','10','location','northeast')
%         
%         

%         figure
%         plot(sum(learns_history))
%         title('# of learning agents')
%         xlabel('time')

        figure
        for k=1:ni
            subplot(1,ni,k)
            area(ind_counts(:,:,k));
            colormap(gca,ind_map)
            title('indicators')
            xlabel('time')
        end
        
        figure
        for k=1:nt
            subplot(1,nt,k)
            area(trait_counts(:,:,k));
            colormap(gca,trait_map)
            title('traits')
            xlabel('time')
        end

        figure
        for k=1:ni
            subplot(1,ni,k)
            plot(ind_shares(k,:));
            colormap(gca,ind_map)
            title('indicators')
            xlabel('time')
        end
        ylim([0 1])
        
        figure
        for k=1:nt
            subplot(1,nt,k)
            plot(trait_shares(k,:));
            colormap(gca,trait_map)
            title('traits')
            xlabel('time')
        end
        ylim([0 1])
        
        for j=1:ni
            for k=1:nt
                figure
                for i=1:range
                   subplot(range,1,i)
                   area(trait_counts_per_inds{i,j}(:,:,k))
                   colormap(trait_map)
                end
                title(sprintf('Traits %d in inds %d',k,j))
            end
        end
        
        
%         figure
%         plot(avg_local_ind_entropy,'r')
%         xlabel('Avg. local ind entropy')

        
        
    end
    
    %toc;
end
