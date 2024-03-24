function [ContStruc,K,best_rho,feas]=Optimize_ConStruc(A,Bdec,Cdec,N)
%Matrix Randomizer 5x5 for Continuous Time keeps on randomizing the matrix and computing gains until rho stops
%decrementing for 100 iterations. Note this is based on the function
%LMI_CT_DeDicont which has to be changed based on the requested control
%algorithm (eg: H2,Hinf,LQ,...)

%first iteration
unique_matrices={};

M = randi([0 1],5,5);   %create random matrix of zero and ones
Diag=diag(ones(5,1));   
ContStruc= M | Diag;%ensure ones on the diagonal

% Convert matrix to a string for comparison
str_matrix = mat2str(ContStruc);
unique_matrices{end+1} = str_matrix;

[K,rho,feas]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStruc);

best_rho = round(rho,5);
counter=0;

while (counter<=100) %change this threshold to reduce/increment iterations after rho stops decreasing
    M_new = randi([0 1],5,5);    
    matrix = M_new | Diag; 
    str_matrix = mat2str(matrix);
        % Check if the generated matrix already exists
    if ~any(strcmp(unique_matrices, str_matrix))
        % If not, add it to the cell array
        unique_matrices{end+1} = str_matrix;
        ContStruc_new = matrix;
    end

    if ContStruc_new == ContStruc
        disp('     ##########   error in generation    #########');
        ContStruc
        ContStruc_new
        break
    else 
        [new_K,new_rho,new_feas]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStruc_new);
        if round(new_rho,5)>=best_rho 
            counter=counter+1
        end
        if round(new_rho,5)<best_rho
            counter=0;
            best_rho=round(new_rho,5);
            feas=new_feas;
            K=new_K;
            ContStruc=ContStruc_new;
        end    
    end
end
