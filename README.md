function f=initialize_variables(N,M,V,min_range,max_range)
  min=min_range;
  max=max_range;
  K=M+V;
  for i=1:N
        for j=1:V
		    f(i,j)=min(j)+(max(j)-min(j))*rand(1);
		end
		f(i,V+1:K)=evaluate_objective(f(i,:),M,V);
  end
end

function f=non_domination_sort_mod(x,M,V)
    [N,~]=size(x);
	clear m
	front=1;
	F(front).f=[];
	individual=[];
	for i=1:N
	    individual(i).n=0;
		individual(i).p=[];
	    for j=1:N
		    dom_less=0;
			dom_equal=0;
			dom_more=0;
			for k=1:M
			    if(x(i,V+k)<x(j,V+k))
				    dom_less=dom_less+1;
				elseif(x(i,V+k)==x(j,V+k))
				    dom_equal=domequal+1;
				else 
				    dom_more=dom_more+1;
				end
			end
			if(dom_less==0&&dom_equal~=M) 
                individual(i).n=individual(i).n+1;
			elseif(dom_more==0&&dom_equal~=M)
			    individual(i).p=[individual(i).p j];
			end
		end
		if(individual(i).n==0)
		    x(i,M+V+1)=1;
			F(front).f=[F(front).f i];
		end
	end
	while ~isempty(F(front).f)
	    Q=[];
		for i=1:length(F(front).f)
		    if(~isempty(individual(F(front).f(i)).p))
			    for j=1:length(individual(F(front).f(i)).p)
		    individual(individual(F(front).f(i)).p(j)).n=individual(individual(F(front).f(i)).p(j)).n-1;					
      if(individual(individual(F(front).f(i)).p(j)).n==0)					         
            x(individual(F(front).f(i)).p(j),M+V+1=front+1;
						Q=[Q individual(F(front).f(i)).p(j)];
					end
				end
			end
		end
		front=front+1;
		F(front).f=Q；
	end
	[temp,index_of_fronts]=sort(x(:,M+V+1));
	for i=1:length(index_of_fronts)
	    sorted_based_on_front(i,:)=x(index_of_front(i),:);
	end
	current_index=0;
	for front=1:(length(F)-1)
	    distance=0;
		y=[];
		previous_index=current_index+1;
		for i=1:length(F(front).f)
		    y(i,:)=sorted_based_on_front(current_index+i,:)
		end
		current_index=current_index+i;%%%current_index=i
		sorted_based_on_objective=[];
		for i=1:M
		    [sorted_based_on_objective,index_of_objectives=sort(y(:,V+i));
			sorted_based_on_objective=[];
			for j=1:length(index_of_objectives)		    
        sorted_based_on_objective(j,:)=y(index_of_objectives(j),:);
				f_max=sorted_based_on_objective(length(index_of_objectives),V+i);
				f_min=sorted_based_on_objective(1,V+i);
				y(index_of_objectives(length(index_of_objectives)),M+V+1+i)=Inf;
				y(index_of_objectives(1),M+V+1+i)=Inf;
				for j=2:length(index_of_objectives)-1 				    
          next_obj=sorted_based_on_objective(j+1,V+i);
					previous_obj=sorted_based_on_objective(j-1,V+i);
					if(f_max-f_min==0)
					    y(index_of_objectives(j),M+V+1+i)=Inf;
					else				    
              y(index_of_objectives(j),M+V+1+i)=(next_obj-previous_obj)/(f_max-f_min);
					end
				end
		end
		distance=[];
		distance(:1)=zeros(length(F(front).f),1);
		for i=1:M
			distance(:,1)=distance(:,1)+y(:,M+V+1+i);
		end
		y(:M+V+2)=distance;
		y=y(:,1:M+V+2);
		z(previous_index:current_index,:)=y;
	end
	f=z();
	
function f=tournament_selection(chromosome,pool_size,tour_size)
    [pop,variables]=size(chromosome);
    rank=variables-1;
    distance=variables;
        for i=1:pool_size  
	    for j=1:tour_size
		    candidate(j)=round(pop*rand(1));
			if(candidate(j)==0)
			    candidate(j)=1;
			end
			if(j>1)
			    while ~isempty(find(candidate(1:j-1)==candidate(j)))
				    candidate(j)=round(pop*rand(1));
					if(candidate(j)==0
					    candidate(j)=1;
					end
				end
			end
		end
		for j=1:tour_size 
		    c_obj_rank(j)=chromosome(candidate(j),rank);
			c_obj_distance(j)=chromosome(candidate(j),distance);
		end
		min_candidate=find(c_obj_rank==min(c_obj_rank));
		if(length(min_candidate)~=1)
		 max_candidate=find(c_obj_distance(min_candidate)==max(c_obj_distance(min_candidate)));
			if(length(max_candidate)~=1
			    max_candidate=max_candidate(1);
			end
			f(i,:)=chromosome(candidate(min_candidate(max_candidate)),:);
		else
	        f(i,:)=chromosome(candidate(min_candidate(1)),:);
		end
	end
end

function f=genetic_operator(parent_chromosome,M,V,mu,mum,l_limit,u_limit)
    [N,m]=size(parent_chromosome);
	clear m
	p=1;
	was_crossover=0;
	was_mutation=0;
	for i=1:N  
	    if(rand(1)<0.9)
		    child_1=[];
			child_2=[];
			parent_1=round(N*rand(1));
			if(parent_1<1)
				parent_1=1;
			end
			parent_2=round(N*rand(1));
			if(parent_2<1)
			    parent_2=1;
			end
			while isequal(parent_chromosome(parent_1,:),parent_chromosome(parent_2,:))
			    parent_2=round(N*rand(1));
				if(parent_2<1)
				    parent_2=1;
				end
			end
			parent_1=parent_chromosome(parent_1,:);
			parent_2=parent_chromosome(parent_2,:);
			for j=1:V
			    u(j)=rand(1);
				if(u(j)<=0.5)
				    bq(j)=(2*u(j))^(1/(mu+1));
				else 
				    bq(j)=(1/(2*(1-u(j))))^(1/(mu+1));
				end
				child_1(j)=0.5*(((1+bq(j))*parent_1(j))+(1-bq(j))*parent_2(j));
				child_2(j)=0.5*(((1-bq(j))*parent_1(j))+(1+bq(j))*parent_2(j));
				if(child_1>u_limit(j))
				    child_1(j)=u_limit(j);
				elseif(child_1(j)<l_limit(j))
				    child_1(j)=l_limit(j);
				end
				if(child_2>u_limit(j))
				    child_2(j)=u_limit(j);
				elseif(child_2(j)<l_limit(j))
				    child_2(j)=l_limit(j);
				end
			end
			child_1(:,V+1:M+V)=evaluate_objective(child_1,M,V);
			child_2(:,V+1:M+V)=evaluate_objective(child_2,M,V);
			was_mutation=0;
			was_crossover=1;
		else
			parent_3=round(N*rand(1));
			if(parent_3<1)
				parent_3=1;
			end
			child_3=parent_chromosome(parent_3,:);
			for j=1:V
			    r(j)=rand(1);
				if(r(j)<0.5)
					delta(j)=(2*r(j))^(1/(mun+1)-1;
				else
					delta(j)=1-(2*(1-r(j)))^(1/(mun+1));
				end
				child_3(j)=child_3(j)+delta(j);
				if(child_3(j)>u_limit(j))
					child_3(j)=u_limit(j);
				elseif(child_3(j)<l_limit(j))
					child_3(j)=l_limit(j);
			    end
			end
			child_3(:,V+1:M+V)=evaluate_objective(child_3,M,V);
			was_mutation=1;
            was_crossover=0;
			if(was_crossover==1)
			    child(p,:)=child_1;
				child(p+1,:)=child_2;
				was_cossover=0;
				p=p+2;
			elseif(was_mutation==1)
			    child(p,:)=child_3(1,1:M+V);
				was_mutation=0;
				p=p+1;
			end
		end
	end
	f=child;
	
function f=replace_chromosome(intermediate_chromosome,M,V,pop)
  [N,m]=size(intermediate_chromosome);
	[temp,index]=sort(intermediate_chromosome(:,M+V+1));
	clear temp m
	for i=1:N
	    sorted_chromosome(i,:)=intermediate_chromosome(index(i),:);
	end
	max_rank=max(intermediate_chromosome(:,M+V+1));
	previous_index=0;
	for i=1:max_rank
	    current_index=max(find(sorted_chromosome(:,M+V+1)==i));
		if(current_index>pop)
		    remaining=pop-previous_index;
			temp_pop=sorted_chromosome(previous_index+1:current_index,:);
			[temp_sort,temp_sort_index]=sort(temp_pop(:,M+V+2),'descend');
			for j=1:remaining
			    f(previous_index+j,:)=temp_pop(temp_sort_index(j),:);
			end
			return 
		elseif(current_index<pop)
		    f(previous_index+1:current_index,:)=sorted_chromosome(previous_index+1:current_index,:);
		end
		previous_index=current_index;
	end
end

function f=evaluate_objective(x,M,V)
    f=[];
	f(1)=x(1);
	g=1;
	sum=0;
	for j=2:V
	    sum=sum+x(i);
	end
	sum=9*(sum/(V-1));
	g=g+sum;
	f(2)=g*(1-sqrt(x(1)/g));
end
