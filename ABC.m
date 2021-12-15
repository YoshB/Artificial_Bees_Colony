% Artificial Bees Colony (ABC)

clc;
clear all
close all

% Choose the function (Ej. 0, 1, 2... etc.)
f_objetivo = 1;

switch f_objetivo
    case 0
        %Griewank
       f = @(x,y) ((x.^2 + y.^2)/4000)-cos(x).*cos(y/sqrt(2)) + 1;
       U = [10 10];
       L = [-10 -10];
       %Minimum = 0; x= 0, y= 0
       
    case 1
        %Rastrigin
       f = @(x,y)  10*2 + x.^2 - 10 .*cos(pi*x)+ y.^2 - 10 .* cos(pi*y);
       U = [5 5];
       L = [-5 -5];
       %Minimum = 0; x = 0, y = 0 
       
      
    case 2
         %DropWave
    f = @(x,y) - ((1 + cos(12*sqrt(x.^2+y.^2))) ./ (0.5 *(x.^2+y.^2) + 2));
    U = [2 2];
    L = [-2 -2];
       %Minimum = -1; x = 0, y =0
       
    case 3
        f = @(x,y) (x-2).^2 + (y-2).^2;
       U = [5 5];
    L = [-5 -5];
    otherwise
        disp("Introduce a valid value of functions")
        return
end

%Iterations
G = 200;
%Population
pob = 50;
%Number of dimensions of the problem, (X,Y) in this case
D = 2;

%Number of attempts
n = 80;

%Worker bees
Pf = 35;
% Observer bees
Po = pob - Pf;

x = zeros(D,Pf);
intentos = zeros(1,Pf);

%initialization
x(1,:) = L(1) + (U(1) - L(1)).*rand(1,Pf);
x(2,:) = L(2) + (U(2) - L(2)).*rand(1,Pf);

%Fitness of the solutions
aptitud = zeros(1,Pf);


[X,Y] = meshgrid(L(1):0.25:U(1), L(2):0.25:U(2));
Z = f(X,Y);
contour(X,Y,Z,35);
hold on

for i=1:G
    
    
    % Workers phase
    for j=1:Pf
        ki = j;
        
        while ki ~= j
            ki = randi([1,Pf]);
        end
        
         % Random between [-1,1]
        phi = -1 + 2*rand;
        
        d = randi([1,D]);
        v = x(:,j);
        % Calculating a new solution
        v(d) = x(d,j) + (x(d,j) - x(d,ki))*phi;
        
        % evaluating if it is better solution
        if f(v(1),v(2)) < f(x(1,j),x(2,j))
            x(:,j) = v;
            intentos(j) = 0;
        else
            %else, let´s sum 1 to its attempts
            intentos(j) = intentos(j) + 1;
        end
    end
    
    %Evaluating every solution and saving their fitness in th aptitude
    %array
    for k=1:Pf
        res = f(x(1,k),x(2,k));

        if res >= 0
            aptitud(k) = 1/(1+res);
        else
            aptitud(k) = 1 + abs(res);
        end
    end
    aptitud = aptitud/sum(aptitud);
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    % Observer Bees Phase
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     

    for j=1: Po
        
        %Choosing a bee based on its aptitude
        m = Ruleta(aptitud);
        
        k = j;
        
        while k ~= j
            k = randi([1,Pf]);
        end
        
        
        d = randi([1,D]);
        phi = -1 + 2*rand;
        
        v = x(:,m);
        %Calculating a candidate solution
        v(d) = x(d,m) + (x(d,m) - x(d,k)) *phi;
        
        %Vif it is a better solution
         if f(v(1),v(2)) < f(x(1,m),x(2,m))
            x(:,m) = v;
            intentos(m) = 0;
         else
            %else, sum 1 to the attempts of the individual
            intentos(m) = intentos(m) + 1;
         end
        
    end
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    % Explorer Phase
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
  
   for j=1:Pf
       %Evaluating if some bee has exceeded its attempts
       %if so, we asign it a random vector
       if intentos(j) > n
           x(1,j) = L(1) + (U(1) - L(1))*rand;
           x(2,j) = L(2) + (U(2) - L(2))*rand;
           intentos(j) = 0;
       end
   end
   
   %Graphics
    axis([L(1),U(1),L(2),U(2)]);
     h = plot(x(1,:), x(2,:),'bx');
     g = plot(x(1,:), x(2,:),'ro');
     pause(0.1)
     delete(h);
     delete(g);
    
end

 axis([L(1),U(1),L(2),U(2)]);
 h = plot(x(1,:), x(2,:),'rx');
 
 aptitud = f(x(1,:), x(2,:));
 
 %Find the minimum of the solutions
[values, y] = min(aptitud)

%The better Solution
x(:,y)
f(x(1,y),x(2,y))



%Function  Ruleta: A function to choose randomly an element based on
%probabilities of each element to be choose
% Input: a list with probabilities of the elements
% Output: the index of the element choosen
function padre= Ruleta(aptitud)

% Guardamos el total de individuos que hay 
 n = length(aptitud);
 padre =1;
     %calculating a random number between 0 and 1
     r = rand;
     p_sum = 0;
     
     for j=1:n
         %sum the probabilitie
        p_sum = p_sum + aptitud(j);
        
        if p_sum >= r
            %selecting this element if the probability is bigger
            padre = j;
            break;;
        end
        
     end
 end
    
