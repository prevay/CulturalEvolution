function [offspring,trait_changes,ind_changes] = create_offspring(range,numagents,ni,nt,np,extra,r,agents,ind_changes,trait_changes,neighborhood,config,n,adjmat)
%Evolutionary algorithm for the cultural evolution model
%  Selects parents and creates offspring agents through recombination and
%  mutation

   %preallocation
   offspring_inds=randi(range,numagents,ni);
   offspring_learns=zeros(numagents,1);
   offspring_traits=randi(range,numagents,nt);
   offspring_fitness=zeros(numagents,1);
   offspring_goodints=zeros(numagents,1);
   offspring_ints=zeros(numagents,1);
   offspring_inds_interacted=zeros(numagents,range,ni);
   if strcmp(config(1),'1')
        offspring_learnfreq=zeros(numagents,1);
   end
   if strcmp(config(2),'1')
    offspring_preft=zeros(numagents,range,np);
    offspring_prefn=zeros(numagents,range,np);
    offspring_badn=zeros(numagents,range,np);
    offspring_badt=zeros(numagents,range,np);
    if strcmp(config(1),'1');
       offspring_learnedn=zeros(numagents,range,nt);
       offspring_learnedt=zeros(numagents,range,nt);
    end
   else
    offspring_prefs=zeros(numagents,range,extra,np);
    if strcmp(config(1),'1');
       offspring_learned_traits=zeros(numagents,extra,nt); 
    end
   end 
   if strcmp(config(4),'1') && strcmp(config(5),'1')
      offspring_p2i=zeros(numagents,2); 
   end
   
   %std_preft=std(agents.pref_t(~isinf(agents.pref_t)));
   %max_preft=max(agents.pref_t(~isinf(agents.pref_t)));
   %std_prefn=mean(std(agents.pref_n));

   
   for i=1:numagents
        %define the neighborhood of each agent:
        if strcmp(config(3),'1') %grid case 
            nbhd=[agents.nbhd(i,:) i];
        elseif strcmp(config(8),'1')
            nbhd=[find(adjmat(i,:)) i]; %network case
        else
            nbhd=1:numagents; %complete case
        end
        
        %pick parents in tournaments of size 3 (from nbhd) (or less if nbhd size <3)
        if length(nbhd)==1 %special case when nbhd size = 1
            parent_idx(1)=nbhd;
            parent_idx(2)=nbhd;
        elseif length(nbhd)==2 %special case for nbhd size = 2
            parent_idx(1)=nbhd(1);
            parent_idx(2)=nbhd(2);
        elseif length(nbhd)==3; %special case for nbhd size = 3
            [~,most_fit]=max(agents.fitness(nbhd));
            parent_idx(1)=nbhd(most_fit);
            nbhd(most_fit)=[];
            [~,most_fit]=max(agents.fitness(nbhd));
            parent_idx(2)=nbhd(most_fit);
        else %all other cases:
            k=1;
            while k<3
                candidates=randsample(nbhd,3);
                if max(candidates)>length(agents.fitness)
                   sprintf('woah') 
                end
                [~,most_fit]=max(agents.fitness(candidates));
                parent_idx(k)=candidates(most_fit);
                while k==2 && parent_idx(1)==parent_idx(2) 
                    k=1;
                end
                k=k+1;
            end
        end

        %RECOMBINATION AND MUTATION:
        
        
        %learning frequency
        if strcmp(config(1),'1')
            if rand<0.01
               offspring_learnfreq(i)=rand;
            elseif rand<0.5
               offspring_learnfreq(i)=agents.learnfreq(parent_idx(1));
            else
               offspring_learnfreq(i)=agents.learnfreq(parent_idx(2)); 
            end
            
            %for ACT-R scenario
            if strcmp(config(2),'1')
                for k=1:range
                    for j=1:nt
                        if rand<0.01
                           offspring_learnedn(i,k,j)=randi(20*n,1);
                        elseif rand<0.5
                           offspring_learnedn(i,k,j)=agents.learned_n(parent_idx(1),k,j);
                        else
                           offspring_learnedn(i,k,j)=agents.learned_n(parent_idx(2),k,j); 
                        end

                        if rand<0.01
                           offspring_learnedt(i,k,j)=randi(20*n*numagents,1);
                        elseif rand<0.5
                           offspring_learnedt(i,k,j)=agents.learned_t(parent_idx(1),k,j);
                        else
                           offspring_learnedt(i,k,j)=agents.learned_t(parent_idx(2),k,j); 
                        end

                    end
                end
            else
                for j=1:nt
                    if rand<0.01
                       offspring_learned_traits(i,j)=randi(range,1);
                    elseif rand<0.5
                       offspring_learned_traits(i,j)=agents.learned_traits(parent_idx(1),j);
                    else
                       offspring_learned_traits(i,j)=agents.learned_traits(parent_idx(2),j); 
                    end 
                end
            end
        end

        %1% mutation and element-wise crossover from two parents:
        %tags:
        for j=1:ni
            if rand<0.01
               offspring_inds(i,j)=randi(range,1);
            elseif rand<0.5
               offspring_inds(i,j)=agents.indicators(parent_idx(1),j);
            else
               offspring_inds(i,j)=agents.indicators(parent_idx(2),j); 
            end
            if strcmp(config(3),'1') && offspring_inds(i,j)~=agents.indicators(i,j);
                ind_changes(i,j)=ind_changes(i,j)+1;
            end
        end
       
        %traits
        for j=1:nt
            if rand<0.01
               offspring_traits(i,j)=randi(range,1);
            elseif rand<0.5
               offspring_traits(i,j)=agents.traits(parent_idx(1),j);
            else
               offspring_traits(i,j)=agents.traits(parent_idx(2),j);
            end
            if strcmp(config(3),'1') && offspring_traits(i,j)~=agents.traits(i,j);
                trait_changes(i,j)=trait_changes(i,j)+1;
            end
        end
        
        %ACT-R preferences:
        if strcmp(config(2),'1')    
            for k=1:range
                 for j=1:np
                            if rand<0.01
                               %if isinf(offspring_preft(i,k,j)) 
                                %  offspring_preft(i,k,j)=max_preft+(-1)^(randi(2,1))*std_preft; 
                               %else
                                %  offspring_preft(i,k,j)=offspring_preft(i,k,j)+(-1)^(randi(2,1))*std_preft; 
                               %end
                               offspring_preft(i,k,j)=randi(20*n*numagents,1);
                            elseif rand<0.5
                               offspring_preft(i,k,j)=agents.pref_t(parent_idx(1),k,j);
                            else
                               offspring_preft(i,k,j)=agents.pref_t(parent_idx(2),k,j); 
                            end

                            if rand<0.01
                               %offspring_prefn(i,k,j)=offspring_prefn(i,k,j)+(-1)^(randi(2,1))*std_prefn; 
                               offspring_prefn(i,k,j)=randi(20*n,1);
                            elseif rand<0.5
                               offspring_prefn(i,k,j)=agents.pref_n(parent_idx(1),k,j);
                            else
                               offspring_prefn(i,k,j)=agents.pref_n(parent_idx(2),k,j); 
                            end

                            if rand<0.01
                               offspring_badt(i,k,j)=randi(20*n*numagents,1);
                            elseif rand<0.5
                               offspring_badt(i,k,j)=agents.bad_t(parent_idx(1),k,j);
                            else
                               offspring_badt(i,k,j)=agents.bad_t(parent_idx(2),k,j); 
                            end

                            if rand<0.01
                               offspring_badn(i,k,j)=randi(20*n,1);
                            elseif rand<0.5
                               offspring_badn(i,k,j)=agents.bad_n(parent_idx(1),k,j);
                            else
                               offspring_badn(i,k,j)=agents.bad_n(parent_idx(2),k,j); 
                            end
                 end
            end
        %NAIVE preferences:    
        else
           for k=1:extra
               for j=1:np
                    if rand<0.01
                       offspring_prefs(i,k,j)=randi(range,1);
                    elseif rand<0.5
                       offspring_prefs(i,k,j)=agents.prefs(parent_idx(1),k,j);
                    else
                       offspring_prefs(i,k,j)=agents.prefs(parent_idx(2),k,j); 
                    end

                    if rand<0.01
                       offspring_preflens(i,k,j)=randi(floor(r*numagents),1);
                    elseif rand<0.5
                       offspring_preflens(i,k,j)=agents.preflens(parent_idx(1),k,j);
                    else
                       offspring_preflens(i,k,j)=agents.preflens(parent_idx(2),k,j); 
                    end
               end
           end
        end
        
        %preference to tags mappings:
        if strcmp(config(4),'1') && strcmp(config(5),'1')
           for k=1:2
               if rand<0.01
                      offspring_p2i(i,k)=randi(2,1);
               elseif rand<0.5
                      offspring_p2i(i,k)=agents.p2i(parent_idx(1),k);
               else
                      offspring_p2i(i,k)=agents.p2i(parent_idx(2),k); 
               end
           end
        end
        
        


   end      
   
   
   %put it all together:
   offspring=struct('indicators',offspring_inds,'traits',offspring_traits,'fitness',...
       offspring_fitness,'good_interactions',offspring_goodints,'interactions',...
       offspring_ints,'indicators_interacted',offspring_inds_interacted);
   if strcmp(config(1),'1')
      offspring.learnfreq=offspring_learnfreq;
      if strcmp(config(2),'1')
         offspring.learned_n=offspring_learnedn;
         offspring.learned_t=offspring_learnedt;
      else
         offspring.learned_traits=offspring_learned_traits;
      end
   end
   
   if strcmp(config(2),'1')
      %offspring.learns=offspring_learns;
      offspring.pref_n=offspring_prefn;
      offspring.pref_t=offspring_preft;
      offspring.bad_n=offspring_badn;
      offspring.bad_t=offspring_badt;

   else
      offspring.prefs=offspring_prefs;
      offspring.preflens=offspring_preflens;
   end
   
   if strcmp(config(3),'1')
      offspring.nbhd=neighborhood; 
   end
   
   if strcmp(config(4),'1') && strcmp(config(5),'1')
      offspring.p2i=offspring_p2i;
   end


end

