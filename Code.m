%A binary-coded GA with Roulette-wheel reproduction scheme,
%two point crossover and bit-wise mutation
%---------------------------------------------------------
clc
clear all

pop=input('Enter population size: ');      %population size
str=40;                                    %string length
pc=input('Enter crossover probability: '); %crossover probability
pm=input('Enter mutation probability: ');  %mutation probability
x1_min=0;                                  %minimum and maximum values of x1 and x2
x1_max=0.5;
x2_min=0;
x2_max=0.5;

itn=0;                                     %iteration
x=randi([0 1],pop,str);                    %creation of first generation
while (itn<=500)                            
for i=1:pop                                %dividing the string length into x1 and x2       
    for j=1:str
        if (j<=str/2)
            x1(i,j)=x(i,j);                %matrix of substring of x1 
        else
            x2(i,j-(str/2))=x(i,j);        %matrix of substring of x2
        end
    end
end

decx1=zeros(pop,1);
decx2=zeros(pop,1);

%decoding the values
for i=1:pop
    for j=1:str/2
        decx1(i,1)=decx1(i,1)+x1(i,j)*2^(str/2-j);
        decx2(i,1)=decx2(i,1)+x2(i,j)*2^(str/2-j);
    end
end 

x1_val=zeros(pop,1);
x2_val=zeros(pop,1);

for i=1:pop
    x1_val(i,1)=x1_val(i,1)+x1_min+(((x1_max-x1_min)/(2^(str/2)-1))*decx1(i,1));
    x2_val(i,1)=x2_val(i,1)+x2_min+(((x2_max-x2_min)/(2^(str/2)-1))*decx2(i,1));
end


x1_real(itn+1)=mean(x1_val(:,1));
x2_real(itn+1)=mean(x2_val(:,1));


%fit=zeros(pop,1);

%for i=1:pop
    %fit(i,1)=fit(i,1)+x1_val(i,1)+x2_val(i,1)-(2*(x1_val(i,1))^2)-(x2_val(i,1)^2)+x1_val(i,1)*x2_val(i,1);
%end


%fitness evaluation
fitx=zeros(pop,1);
for i=1:pop
    fitx(i,1)=func_eval(x1_val(i,1),x2_val(i,1));
end

min_fitness(itn+1)=min(fitx);
max_fitness(itn+1)=max(fitx);
avg_fitness(itn+1)=sum(fitx)/pop;

fitnew=fitx/sum(fitx);

%selection of mating pool using Roulette Wheel Selection
for i=2:pop
    fitnew(i)=fitnew(i)+fitnew(i-1);
end

for i=1:pop
    random=rand;
    for j=1:pop
        if random<=fitnew(j)
            matingpool(i,:)=x(j,:);
            break
        end
    end
end

    randomx=randperm(pop,pop);
    for i=1:length(randomx)
        newmat(i,:)=matingpool(randomx(i),:);
    end
    
%two-point crossover    
    for i=1:2:pop
        randomr=rand;
    if randomr<pc
        gen=randperm(str/2-1,2);
        if gen(1)<gen(2)
            store=newmat(i,gen(1):gen(2));
            newmat(i,gen(1):gen(2))=newmat(i+1,gen(1):gen(2));
            newmat(i+1,gen(1):gen(2))=store;
            else
            store=newmat(i,gen(2):gen(1));
            newmat(i,gen(2):gen(1))=newmat(i+1,gen(2):gen(1));
            newmat(i+1,gen(2):gen(1))=store;
        end 
    end
    end

%mutation
for i=1:pop
    for j=1:str
        randomz=rand;
        if randomz<pm
            if newmat(i,j)==0
               newmat(i,j)=1;
            else
               newmat(i,j) = 0;
            end
        end
    end
end

x=newmat;
itn=itn+1;

end

%plotting
figure(1)
plot(1:itn,avg_fitness);
axis([1 itn 0.5 1]);
xlabel('Generations');
ylabel('Average fitness');
legend('Average fitness');
title('Average fitness vs No. of generations');
hold off;


figure(2)
plot(1:itn,max_fitness);
axis([1 itn 0.5 1]);
xlabel('Generations');
ylabel('Fitness');
hold on;

plot(1:itn,min_fitness);
axis([1 itn 0.5 1]);
xlabel('Generations');
ylabel('Fitness');
legend('Maximum fitness','Minimum fitness');
title('Minimum fitness and Maximum fitness vs No. of generations');
hold off;

figure(3)
plot(1:itn,x1_real);
axis([1 itn 0 1]);
xlabel('Generations');
ylabel('Variable values');
hold on;

plot(1:itn,x2_real);
axis([1 itn 0 1]);
xlabel('Generations');
ylabel('Variable values');
legend('x1','x2');
title('Optimal Solution');
hold off;


%md = max(fitx);
%minimumfunc = ((1/md)-1)^0.5

minfunc_value=(1/max(fitx)-1)^0.5
x1_value=min(x1_val)
x2_value=min(x2_val)


%objective function
function f=func_eval(x1,x2)

eqn=fileread('Input.txt');
fh=str2func(eqn);
f1=fh(x1,x2);
f=1/(1+f1^2);
end