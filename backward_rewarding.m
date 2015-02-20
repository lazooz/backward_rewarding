% This script was written for Lazooz to have a backward rewarding
% algorithm.
clear all


% % % Reading the matrix
[num1,txt1,raw1] = xlsread('Data_delete2.xlsx'); % enter the excel file
num1_no_nan = num1;
num1_no_nan(isnan(num1)) = 0;
num1_temp = num1;
txt_replaced =cell(length(num1),1);
txt_replaced(1:length(num1))= txt1(2:length(num1)+1); 
txt_replaced_temp =cell(length(num1),1);
txt_replaced_temp(1:length(num1))= txt1(2:length(num1)+1); 
% % % Filling the slots of people who reccomended other people


for yyy = 1: length(num1)
for ttt = 1: length(num1)
    
    for jjj = 1:length(num1)      
    if isequal(txt1{yyy+1,ttt+1},txt1{1,jjj+1}) == 1 && isnan(num1(jjj,ttt)) == 0 && isnumeric(num1(jjj,ttt)) == 1 %This is the case someone reccomended on someone else. checking  that the recomended person evaluated
   
        mutual_index = max(abs(isnan(num1(yyy,:))-1) + abs(isnan(num1(jjj,:))-1) -1,0);  %1 if both people estimated a person so there is a reference
         tot_mutual_gave = num1_no_nan(yyy,:) *mutual_index';%the person gave his right to the other
        tot_mutual_eval = num1_no_nan(jjj,:) *mutual_index'; % the person who has the authoroty
        
num1_temp(yyy,ttt)  = num1(jjj,ttt).*tot_mutual_gave./tot_mutual_eval;

    end
    end
  
end
end
    
    num1_temp_no_nan = num1_temp;
    num1_temp_no_nan(isnan(num1_temp)) = 0; 
rescale_factor = repmat(sum(num1_temp_no_nan,2),1,length(num1)); 
num1_temp_no_nan_rescaled = num1_temp_no_nan./rescale_factor;

num1_temp_rescaled = num1_temp_no_nan_rescaled;
num1_temp_rescaled(num1_temp_rescaled==0) = nan;




A = num1_temp_rescaled;

C = A; 


% I want to make sure it doesn't matter to change the order
%This is also a scheme how to switch between rows. in case I have a problematic order in the beginning 
%%%%%%%%%%%%%%%%%%%%%%%%
% { 
for ooo = 1:length(num1)

    while mean(isnan(A(1:ooo,ooo))) ==  1 || mean(isnan(A(ooo,1:ooo-1))) ==  1   %check if I need to reorder if there is no common evluation
    change_index1 = ooo; %choosing which raws to switch. I do it in a random process here so I won't be stuck in repititions
change_index2 = max(ceil(rand(1)*length(A)),ooo+1); 

C = A;
A(change_index1,:) = C(change_index2,:);
A(:,change_index1) = C(:,change_index2);

A(change_index2,:) = C(change_index1,:);
A(:,change_index2) = C(:,change_index1);

A(change_index2,change_index2)  = C(change_index1,change_index1) ;
A(change_index1,change_index1)  = C(change_index2,change_index2) ;

A(change_index1,change_index2)  = C(change_index2,change_index1) ;
A(change_index2,change_index1)  = C(change_index1,change_index2) ;
txt_replaced_temp(ooo) = txt_replaced(ooo);
txt_replaced(ooo) = txt_replaced(change_index2); %changing the name ordering
txt_replaced(change_index2) = txt_replaced_temp(ooo);
    end    
end

%}
txt_replaced
%%

nan_matrix = isnan(A);
max_index = 1; % need to think how to calculate it when I rearange the matrix
iterations_num = 30 ; % How many iterations to do. From what I played with it seems that 10-20 are enough

for  i = max_index: length(A)
temp_matrix = A(1:i,1:i);
temp_matrix_no_nan = temp_matrix;
temp_matrix_no_nan(isnan(temp_matrix)) = 0;
rescale_factor = repmat(sum(temp_matrix_no_nan,2),1,i); 
temp_matrix_rescaled = temp_matrix_no_nan./rescale_factor;
temp_matrix_rescaled_temp = temp_matrix_no_nan./rescale_factor;



for xxx = 1: iterations_num

if i == max_index %If it is the first itteration I should skip the if since I don't have Eigenvectors yet.
    
else %filling the missing entries
    
        if length(norm_V) == i %check if it is first interaction after added a new person
norm_V_temp =norm_V';

        else    
         norm_V_temp =[norm_V' 0]; %So I could estimate the weight of the new person in the first iteration - at the moment I put zero
        end



    for jjj = 1:i % Check who didn't estimate the new person 
           for lll = 1:i %checking who did you estimate  
 
        %Here I need to make loop since I need to repeat all the
        %people estimation in every iteration
        
if isnan(temp_matrix(jjj,lll )) == 1 % filling the missing entries
    
 temp_matrix_rescaled(jjj,lll ) = norm_V_temp(lll); %entering the value of the eigen vector for the unknown values

end
    
           end
    
        
   
    end
 
    
    % Here I need to normalize the matrix again. I need to normalize only
    % the values that I didn't enter by force. 
         normalize_factor_vector = (sum(temp_matrix_rescaled.*(abs(isnan(temp_matrix)-1)),2))./(1 - sum(temp_matrix_rescaled(:,:).* isnan(temp_matrix(:,:)),2)); %%

         
         
   
         
         % The normalization process doesn't have 1 in the numerator since
         % in every iteration I do the numbers in the original matrix might be
         % different than unity
         normalize_factor_vector_mat = repmat(normalize_factor_vector,1,6);
         not_nan_logical = logical(abs(isnan(temp_matrix(:,:)) - 1));
         temp_matrix_rescaled(not_nan_logical) = temp_matrix_rescaled(not_nan_logical)./...
             normalize_factor_vector_mat(not_nan_logical) ; 
         
         
         
     
    
end



[V,D] = eig(temp_matrix_rescaled'); % check later how I get the data and see that it makes sense if I need to take a prime
max_EV = max(max(D));%should be unity
if max_EV >(1+10^(-6))  || max_EV <(1-10^(-6)) % making sure my eigenvalue is one (otherwise there is a mistake).
disp('max_EV is not equal to 1');
max_EV
end

[row,col] = find(D == max_EV); % making sure I am taking the eigenvector that is related to the eigenvalue
norm_V = abs(V(:,row))./(sum(abs(V(:,row)),1));


end

end

norm_V  %the weights (to this we will have to add the old weights such that the distribution of zooz will take into account both - can discuss with yani if you have questions)
txt_replaced %the name each number represent
