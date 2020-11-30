% nsub is function that returns a non singular submatrix
%
% [R,C] = nsub(A,r) returns r rows and columns that are linear independent
%
% [R,C,time] = nsub(A,r) returns r rows and columns that are linear
% independent with the elapsed time to run the algorithm
%
% A(R,C) is a non-singular submatrix of A
%
% Example:
%   r = 3;
%   rc = linspace(1,10,r);
%   A = full(sprand(5,5,1,rc));
%   [R,C] = nsub(A,r);
%   rank(A(R,C)) == r % Check if it is true

function [R,C,time] = nsub(A,r)

    tic % start timer

    [m,n] = size(A); % get number of rows and cols 
    
    rankInit = rank(A,10^-6);
    
    if rankInit < r
        fprintf(['\nThe instance A does not have ', num2str(r),' singular values well defined, please set r to be less or equal than ', num2str(rankInit),'\n\n']);
        R = [];
        C = [];
        time = toc;
    else
        
    alpha = 1.25; % default

    maxIter = ceil(max(5,alpha*10^(-4)*n*r)); % maximum number of random sets of rows/cols (default)
    
    penValue_Min = 0; % value for the initial penalized identity submatrix (default).
    
    penValue_Max = 4; % value for the final penalized values of the identity submatrix (default).
    
    %% 
    
    if m > n
                
        counter = 0; % counter of random sets of cols

        flag_Row = 0; % to start while cycle

        C_tot = [1:n];  % All columns

        best_sizeR = -1; % Best size of linear independent rows

        while flag_Row == 0 % start while cycle

            counter = counter + 1; % update counter

            if counter > 1 % if counter > 1 -> reorder

                rng('shuffle','twister'); % control random number generator

                C_tot = (randperm(length(C_tot))); % reorder columns

            end

            % Penalized Identity matrix above A
            
            A_new = [zeros(r,n); A]; % new matrix A with the penalized identity submatrix above 
            
            penValue = penValue_Min;
            
            m1 = m + r; % Number of rows of the matrix A with the penalized identity submatrix above

            R  = [1:r]; % initial rows of the submatrix
            Rb = [r+1:m1]; % initial rows that aren't in the submatrix

            C = C_tot(1:r); % columns of the submatrix

            for i = (1:r)
    
                id_R = R(i); % submatrix Row indice
                id_C = C(i); % submatrix Col indice

                A_new(id_R,id_C) = 1* 10^(-penValue); % replace 0 to penalized value to create an Linear Independent submatrix

            end
            
            countRi = 0; % counter of Linear Independent Rows founded

            flag_LSFI = 1; % to start the local search changing rows
            
            while flag_LSFI == 1 && countRi < r
                               
                flag_LSFI = 0;

                Ar = A_new(R,C);

                [L1,U1,P1] = lu(Ar');

                i = 1;

                while i <= m1-r-countRi

                    % LU Factoration   

                    b = P1 * A_new(Rb(i),C)';

                    y = L1\b;

                    alfa = U1\y;

                    % Changing Ar

                    local_alfa = find(abs(alfa)>1);

                    T = isempty(local_alfa);

                    if T == 0

                        [tamAlfa,~] = size(local_alfa);

                        flagAlfa = 0;

                        for iAlfa = (1:tamAlfa)

                            valueAlfa = local_alfa(iAlfa);

                            if R(valueAlfa) <= r
                                
                                local_alfa = valueAlfa; 
                                flagAlfa = 1; 

                                break
                            end
                        end

                        if flagAlfa == 0
                            local_alfa = local_alfa(1);

                        end

                        if R(local_alfa) <= r % Preference to change rows of Penalized Identity submatrix 
                     
                            R(local_alfa) = Rb(i);
                            Rb(i) = []; %Remove element from matrix

                            countRi  = countRi + 1;

                            flag_LSFI = 1;
                            
                            if countRi==r
                                break
                            end

                        else

                            local_alfa = local_alfa(1);

                            el_save = R(local_alfa) ;

                            R(local_alfa) = Rb(i);

                            Rb(i) = el_save;

                        end

                        Ar = A_new(R,C);

                        [L1,U1,P1] = lu(Ar');
                    end
                    i = i + 1;
                end
                
                if flag_LSFI == 0 && countRi < r && penValue < penValue_Max
                    
                    penValue = penValue + 1;
                    
                    for i = (1:r)

                        id_R = R(i); % submatrix Row indice
                        
                        if id_R <= r
                            
                            id_C = C(i); % submatrix Col indice
                            A_new(id_R,id_C) = 1* 10^(-penValue); %update penalized value to create an Linear Independent submatrix
                            
                        end
                    end
                    
                    flag_LSFI = 1;
                    
                end

            end % End of local search
            
            R_init = R; % copy R_init from R
            
            R = []; % Actual linear independent rows of A
            count_R = 0; % number of linear independent rows of A

            for j = (1:r)

                Rj_temp = R_init(j) - r; % A_new has r more elements than A. If this indice is negative or 0, so it doesn't belongs to A

                if Rj_temp > 0 % if it belongs to A

                    R = [R;Rj_temp]; % Add this row to the set of Linear Independent Rows
                    
                    count_R = count_R + 1; % Update linear independent rows counter
                end
            end
            
            if count_R == r % Can finish the algorithm, because it was found the r linear independent rows
                
                flag_Row = 1; % Update flag
                
            end

            if count_R > best_sizeR % if we found more linear independent rows than before -> update the best linear independent rows

                best_sizeR = count_R; % update the maximum number of linear independent rows found

                best_R = R; % best set of linear independent rows
                best_C = C; % best set of linear independent columns

            end
            
            
            if counter == maxIter && flag_Row == 0 % If it reaches the maximum number of iteration and we didn't find r linear independent rows
                
                flag_Row = 2; % Update flag
                
            end
              
        end
        
        R = best_R; % Update Rows
        C = best_C; % Update Cols

        if flag_Row == 2 % If we didn't find r linear independent rows

            minSigmaFactor = abs(floor(log10(abs(min(svd(A(R,:),'econ')))))); % find min factor of singular values of A

            C = (1:n); % Use all columns
            
            % Start Greedy Light 
            
            sizeR = length(R);

            while sizeR < r

                flag = 0;

                sizeRows = length(R);
                sizeRows = sizeRows+1;

                for i = (1:m)

                    X = isempty(find(R==i));

                    if X==1

                        R_save = [R;i];

                        min_sigma = min(svd(A(R_save,:),'econ'));

                        if min_sigma > 10^(-minSigmaFactor) 

                            flag = 1;

                            R = [R;i];
                            break
                        end
                    end
                end

                if flag == 0
                   minSigmaFactor = minSigmaFactor + 1;
                end

                sizeR = length(R);

            end % Finish Greedy Light
            
            % Penalized Identity matrix on the left of A
            
            A_new = [zeros(m,r) A];  % new matrix A with the penalized identity submatrix on the left
            
            penValue = penValue_Min;
            
            for i = (1:r)
    
                id_R = R(i); % submatrix Row indice
                id_C = C(i); % submatrix Col indice

                A_new(id_R,id_C) = 1* 10^(-penValue); % replace 0 to penalized value to create an Linear Independent submatrix

            end
            
            n1 = n + r; % number of columns of the matrix A with the penalized identity submatrix above

            C  = [1:r]; % initial columns of the submatrix

            Cb = [r+1:n1]; % initial rows that aren't in the submatrix
            
            countCi = 0; % counter of Linear Independent Rows founded
            
            flag_Col = 1; % flag for while cycle
            
            % Start local search to change columns
            
            while flag_Col == 1 && countCi < r

                flag_Col = 0;

                Ar = A_new(R,C);

                % FOR COLUMNS

                [L2,U2,P2] = lu(Ar);

                i = 1;

                while i <= n1-r-countCi

                    b = P2 * A_new( R,Cb(i));

                    y = L2\b;

                    alfa = U2\y;

                    % Changing Ar

                    local_alfa = find(abs(alfa)>1);

                    T = isempty(local_alfa);

                    if T == 0

                        [tamAlfa,~] = size(local_alfa);

                        flagAlfa = 0;

                        for iAlfa = (1:tamAlfa)

                            valueAlfa = local_alfa(iAlfa);

                            if C(valueAlfa) <= r %Alteracao v2, nao tinha C no valueAlfa

                                local_alfa = valueAlfa;
                                flagAlfa = 1;

                                break
                            end
                        end

                        if flagAlfa == 0
                            local_alfa = local_alfa(1);

                        end

                        if C(local_alfa) <= r

                            countCi  = countCi + 1;

                            C(local_alfa) = Cb(i);
                            Cb(i) = []; %Remove element from matrix

                            flag_Col = 1;

                            if countCi==r
                                break
                            end

                        else

                            local_alfa = local_alfa(1);

                            el_save = C(local_alfa) ;

                            C(local_alfa) = Cb(i);

                            Cb(i) = el_save;
                        end

                        Ar = A_new(R,C);

                        [L2,U2,P2] = lu(Ar);
                    end
                    i = i + 1;
                end
                
                if flag_Col == 0 && countCi < r && penValue < penValue_Max
                    
                    penValue = penValue + 1;
                    
                    for i = (1:r)

                        id_C = C(i); % submatrix Row indice
                        
                        if id_C <= r
                            
                            id_R = R(i); % submatrix Col indice
                            A_new(id_R,id_C) = 1* 10^(-penValue); %update penalized value to create an Linear Independent submatrix
                            
                        end
                    end
                    
                    flag_Col = 1;
                    
                end

            end % end local search
            
            C_init = C; % copy C_init from C
            
            C = []; % Actual linear independent columns of A
            count_C = 0; % number of linear independent rows of A

            for j = (1:r)

                Cj_temp = C_init(j) - r; % A_new has r more elements than A. If this indice is negative or 0, so it doesn't belongs to A

                if Cj_temp > 0 % if it belongs to A

                    C = [C;Cj_temp]; % Add this row to the set of Linear Independent Rows
                    
                    count_C = count_C + 1; % Update linear independent rows counter
                end
            end
            
            if count_C < r % if it wasn't found all r cols linear independent, flag_Col = 2
                flag_Col = 2;
            end
    
            
            if flag_Col == 2 %Go to Greedy Light if we don't have r linear independent cols 

                minSigmaFactor = abs(floor(log10(abs(min(svd(A(R,C),'econ')))))); % find min factor of singular values of A

                % Start Greedy Light 

                sizeC = length(C);

                while sizeC < r

                    flag = 0;

                    sizeColumns = length(C);
                    sizeColumns = sizeColumns+1;

                    for i = (1:n)

                        X = isempty(find(C==i));

                        if X==1
 
                            C_save = [C;i];

                            min_sigma = min(svd(A(R,C_save),'econ'));

                            if min_sigma > 10^(-minSigmaFactor)

                                flag = 1;

                                C = [C;i];

                                break
                            end
                        end

                    end

                    if flag == 0
                       minSigmaFactor = minSigmaFactor + 1; 
                    end

                    sizeC = length(C);

                end % end Greedy Light
            end % end if didn't find all r linear independent cols
        end % end if didn't find all r linear independent rows
       
        
    elseif m < n
        %% start the case if n > m

        counter = 0; % counter of random sets of rows

        flag_Col = 0; % to start while cycle

        R_tot = [1:m];  % All rows

        best_sizeC = -1; % Best size of linear independent cols

        while flag_Col == 0
            
            counter = counter + 1; % update counter

            if counter > 1 % if counter > 1 -> reorder

                rng('shuffle','twister'); % control random number generator

                R_tot = (randperm(length(R_tot))); % reorder rows

            end
            
            n1 = n + r; % Number of cols s of the matrix A with the penalized identity submatrix above
            
            C  = [1:r]; % initial cols of the submatrix
            
            Cb = [r+1:n1]; % initial cols that aren't in the submatrix

            R = R_tot(1:r); % rows of the submatrix

            % Penalized Identity matrix on the left of A
            
            A_new = [zeros(m,r) A]; % new matrix A with the penalized identity submatrix on the left
           
            penValue = penValue_Min;
            
            for i = (1:r)
    
                id_R = R(i); % submatrix Row indice
                id_C = C(i); % submatrix Col indice

                A_new(id_R,id_C) = 1* 10^(-penValue); % replace 0 to penalized value to create an Linear Independent submatrix

            end

            countCi = 0; % counter of Linear Independent Columns founded

            flag_LSFI = 1; % to start the local search changing columns
            
            while flag_LSFI == 1 && countCi < r % Start local search

                flag_LSFI = 0;
                
                Ar = A_new(R,C);

                % FOR COLUMNS

                [L2,U2,P2] = lu(Ar);

                i = 1;

                while i <= n1-r-countCi

                    b = P2 * A_new( R,Cb(i));

                    y = L2\b;

                    alfa = U2\y;

                    % Changing Ar

                    local_alfa = find(abs(alfa)>1);

                    T = isempty(local_alfa);

                    if T == 0

                        [tamAlfa,~] = size(local_alfa);

                        flagAlfa = 0;

                        for iAlfa = (1:tamAlfa)

                            valueAlfa = local_alfa(iAlfa);

                            if C(valueAlfa) <= r %Alteracao v2, nao tinha C no valueAlfa

                                local_alfa = valueAlfa;
                                flagAlfa = 1;

                                break
                            end
                        end

                        if flagAlfa == 0
                            local_alfa = local_alfa(1);

                        end

                        if C(local_alfa) <= r
                            
                            countCi  = countCi + 1;

                            C(local_alfa) = Cb(i);
                            Cb(i) = []; %Remove element from matrix

                            flag_LSFI = 1;

                            if countCi==r
                                break
                            end
                                
                        else

                            local_alfa = local_alfa(1);

                            el_save = C(local_alfa) ;

                            C(local_alfa) = Cb(i);

                            Cb(i) = el_save;
                        end

                        Ar = A_new(R,C);

                        [L2,U2,P2] = lu(Ar);
                    end
                    i = i + 1;
                end
                
                                
                if flag_LSFI == 0 && countCi < r && penValue < penValue_Max
                    
                    penValue = penValue + 1;
                    
                    for i = (1:r)

                        id_C = C(i); % submatrix Row indice
                        
                        if id_C <= r
                            
                            id_R = R(i); % submatrix Col indice
                            A_new(id_R,id_C) = 1* 10^(-penValue); %update penalized value to create an Linear Independent submatrix
                            
                        end
                    end
                    
                    flag_LSFI = 1;
                    
                end

            end % end local search
            
            C_init = C; % copy C_init from C
            
            C = []; % Actual linear independent cols of A
            count_C = 0; % number of linear independent cols of A

            for j = (1:r)

                Cj_temp = C_init(j) - r; % A_new has r more elements than A. If this indice is negative or 0, so it doesn't belongs to A

                if Cj_temp > 0 % if it belongs to A

                    C = [C;Cj_temp]; % Add this col to the set of Linear Independent Cols
                    count_C = count_C + 1; % Update linear independent cols counter
                end

            end
            
            if count_C == r % Can finish the algorithm, because it was found the r linear independent cols
                
                flag_Col = 1; % Update flag
                
            end

            if count_C > best_sizeC % if we found more linear independent cols than before -> update the best linear independent cols

                best_sizeC = count_C; % update the maximum number of linear independent cols found
                
                best_R = R; % best set of linear independent rows
                best_C = C; % best set of linear independent columns
                
            end

            if counter == maxIter && flag_Col == 0 % If it reaches the maximum number of iteration and we didn't find r linear independent cols
                
                flag_Col = 2; % Update flag
                
            end
            
        end % end flag_Col
 
        R = best_R; % Update Rows
        C = best_C; % Update Cols

        if flag_Col == 2 % If we didn't find r linear independent cols
            
            minSigmaFactor = abs(floor(log10(abs(min(svd(A(:,C),'econ')))))); % find min factor of singular values of A

            R = (1:m); % Use all rows
            
            % Start Greedy Light 

            sizeC = length(C);

            while sizeC < r

                flag = 0;

                sizeColumns = length(C);
                sizeColumns = sizeColumns+1;

                for i = (1:n)

                    X = isempty(find(C==i));

                    if X==1
           
                        C_save = [C;i];

                        min_sigma = min(svd(A(R,C_save),'econ'));

                        if min_sigma > 10^(-minSigmaFactor)

                            flag = 1;

                            C = [C;i];

                            break
                        end
                    end
                end

                if flag == 0
                   minSigmaFactor = minSigmaFactor + 1; 
                end

                sizeC = length(C);

            end % end Greedy Light

            % Penalized Identity matrix above A
            
            A_new = [zeros(r,n); A]; % new matrix A with the penalized identity submatrix above
            
            penValue = penValue_Min;
            
            for i = (1:r)
    
                id_R = R(i); % submatrix Row indice
                id_C = C(i); % submatrix Col indice

                A_new(id_R,id_C) = 1* 10^(-penValue); % replace 0 to penalized value to create an Linear Independent submatrix

            end

            m1 = m + r; % Number of rows of the matrix A with the penalized identity submatrix above

            R  = [1:r]; % initial rows of the submatrix
            Rb = [r+1:m1]; % initial rows that aren't in the submatrix
            
            best_sizeR = -1;

            countRi = 0; % counter of Linear Independent Rows founded

            flag_LSFI = 1; % to start the local search changing rows
            
            while flag_LSFI == 1 && countRi < r

                flag_LSFI = 0;

                Ar = A_new(R,C);

                [L1,U1,P1] = lu(Ar');

                i = 1;

                while i <= m1-r-countRi

                    % LU Factoration   

                    b = P1 * A_new(Rb(i),C)';

                    y = L1\b;

                    alfa = U1\y;

                    % Changing Ar

                    local_alfa = find(abs(alfa)>1);

                    T = isempty(local_alfa);

                    if T == 0

                        [tamAlfa,~] = size(local_alfa);

                        flagAlfa = 0;

                        for iAlfa = (1:tamAlfa)

                            valueAlfa = local_alfa(iAlfa);

                            if R(valueAlfa) <= r %Alteracao v2, nao tinha R no valueAlfa

                                local_alfa = valueAlfa;
                                flagAlfa = 1;

                                break
                            end
                        end

                        if flagAlfa == 0
                            local_alfa = local_alfa(1);

                        end

                        if R(local_alfa) <= r % Preference to change rows of Penalized Identity submatrix 
                            
                            R_save = R;
                            R_save(local_alfa) = Rb(i);
                            
                            if min(svd(A_new(R_save,C),'econ')) > 10^-8 %Alteracao v2

                                countRi  = countRi + 1;

                                R(local_alfa) = Rb(i);
                                Rb(i) = []; %Remove element from matrix

                                flag_LSFI = 1;

                                if countRi==r
                                    break
                                end
                            end

                        else

                            local_alfa = local_alfa(1);

                            el_save = R(local_alfa) ;

                            R(local_alfa) = Rb(i);

                            Rb(i) = el_save;

                        end

                        Ar = A_new(R,C);

                        [L1,U1,P1] = lu(Ar');
                    end
                    i = i + 1;
                end
                                
                if flag_LSFI == 0 && countRi < r && penValue < penValue_Max
                    
                    penValue = penValue + 1;
                    
                    for i = (1:r)

                        id_R = R(i); % submatrix Row indice
                        
                        if id_R <= r
                            
                            id_C = C(i); % submatrix Col indice
                            A_new(id_R,id_C) = 1* 10^(-penValue); %update penalized value to create an Linear Independent submatrix
                            
                        end
                    end
                    
                    flag_LSFI = 1;
                    
                end

            end % End of local search
            
            R_init = R; % copy R_init from R
            
            R = []; % Actual linear independent rows of A
            count_R = 0; % number of linear independent rows of A

            for j = (1:r)

                Rj_temp = R_init(j) - r; % A_new has r more elements than A. If this indice is negative or 0, so it doesn't belongs to A

                if Rj_temp > 0 % if it belongs to A

                    R = [R;Rj_temp]; % Add this row to the set of Linear Independent Rows
                    
                    count_R = count_R + 1; % Update linear independent rows counter
                end
            end
            
            if count_R == r % Can finish the algorithm, because it was found the r linear independent rows
                
                flag_Row = 1; % Update flag
                
            end

            if count_R > best_sizeR % if we found more linear independent rows than before -> update the best linear independent rows

                best_sizeR = count_R; % update the maximum number of linear independent rows found

                best_R = R; % best set of linear independent rows
                best_C = C; % best set of linear independent columns

            end
            
            
            if count_R < r % If we didn't find r linear independent rows
                
                flag_Row = 2; % Update flag
                
            end
              
            R = best_R; % Update Rows
            C = best_C; % Update Cols

            if flag_Row == 2 % If we didn't find r linear independent rows
            
                minSigmaFactor = abs(floor(log10(abs(min(svd(A(R,C),'econ')))))); % find min factor of singular values of A

                % Start Greedy Light 

                sizeR = length(R);

                while sizeR < r

                    flag = 0;

                    sizeRows = length(R);
                    sizeRows = sizeRows+1;

                    for i = (1:m)

                        X = isempty(find(R==i));

                        if X==1
                            
                            R_save = [R;i];

                            min_sigma = min(svd(A(R_save,C),'econ'));

                            if min_sigma > 10^(-minSigmaFactor) 

                                flag = 1;

                                R = [R;i];
                                break
                            end
                        end
                    end

                    if flag == 0
                       minSigmaFactor = minSigmaFactor + 1; 
                    end

                    sizeR = length(R);

                end % end of Greedy Light
            end % end if we didn't find r linear independent rows
        end % end if we didn't find r linear independent cols
      
    else % if m =n
        
        if (sum(sum(abs(A - A'))))<=10^-10 % symmetric case
                    
            counter = 0; % counter of random sets of cols

            flag_Row = 0; % to start while cycle

            C_tot = [1:n];  % All columns

            best_sizeR = -1; % Best size of linear independent rows

            while flag_Row == 0 % start while cycle

                counter = counter + 1; % update counter

                if counter > 1 % if counter > 1 -> reorder

                    rng('shuffle','twister'); % control random number generator

                    C_tot = (randperm(length(C_tot))); % reorder columns

                end

                % Penalized Identity matrix above A

                A_new = [zeros(r,n); A]; % new matrix A with the penalized identity submatrix above 
                
                penValue = penValue_Min;

                m1 = m + r; % Number of rows of the matrix A with the penalized identity submatrix above

                R  = [1:r]; % initial rows of the submatrix
                Rb = [r+1:m1]; % initial rows that aren't in the submatrix

                C = C_tot(1:r); % columns of the submatrix

                for i = (1:r)

                    id_R = R(i); % submatrix Row indice
                    id_C = C(i); % submatrix Col indice

                    A_new(id_R,id_C) = 1* 10^(-penValue); % replace 0 to penalized value to create an Linear Independent submatrix

                end

                countRi = 0; % counter of Linear Independent Rows founded

                flag_LSFI = 1; % to start the local search changing rows

                while flag_LSFI == 1 && countRi < r

                    flag_LSFI = 0;

                    Ar = A_new(R,C);

                    [L1,U1,P1] = lu(Ar');

                    i = 1;

                    while i <= m1-r-countRi

                        % LU Factoration   
                        
                        b = P1 * A_new(Rb(i),C)';

                        y = L1\b;

                        alfa = U1\y;

                        % Changing Ar

                        local_alfa = find(abs(alfa)>1);

                        T = isempty(local_alfa);

                        if T == 0

                            [tamAlfa,~] = size(local_alfa);

                            flagAlfa = 0;

                            for iAlfa = (1:tamAlfa)

                                valueAlfa = local_alfa(iAlfa);

                                if R(valueAlfa) <= r %Alteracao v2, nao tinha R no valueAlfa

                                    local_alfa = valueAlfa;
                                    flagAlfa = 1;

                                    break
                                end
                            end

                            if flagAlfa == 0
                                local_alfa = local_alfa(1);

                            end

                            if R(local_alfa) <= r % Preference to change rows of Penalized Identity submatrix 
                                
                                countRi = countRi + 1;
                                
                                R(local_alfa) = Rb(i);
                                Rb(i) = []; %Remove element from matrix

                                flag_LSFI = 1;

                                if countRi==r
                                    break
                                end
                                
                            else

                                local_alfa = local_alfa(1);

                                el_save = R(local_alfa) ;

                                R(local_alfa) = Rb(i);

                                Rb(i) = el_save;

                            end

                            Ar = A_new(R,C);

                            [L1,U1,P1] = lu(Ar');
                        end
                        i = i + 1;
                    end % end while
                                    
                    if flag_LSFI == 0 && countRi < r && penValue < penValue_Max

                        penValue = penValue + 1;

                        for i = (1:r)

                            id_R = R(i); % submatrix Row indice

                            if id_R <= r

                                id_C = C(i); % submatrix Col indice
                                A_new(id_R,id_C) = 1* 10^(-penValue); %update penalized value to create an Linear Independent submatrix

                            end
                        end

                        flag_LSFI = 1;

                    end
                end % End of local search

                R_init = R; % copy R_init from R

                R = []; % Actual linear independent rows of A
                count_R = 0; % number of linear independent rows of A

                for j = (1:r)

                    Rj_temp = R_init(j) - r; % A_new has r more elements than A. If this indice is negative or 0, so it doesn't belongs to A

                    if Rj_temp > 0 % if it belongs to A

                        R = [R;Rj_temp]; % Add this row to the set of Linear Independent Rows

                        count_R = count_R + 1; % Update linear independent rows counter
                    end
                end

                if count_R == r % Can finish the algorithm, because it was found the r linear independent rows

                    flag_Row = 1; % Update flag

                end

                if count_R > best_sizeR % if we found more linear independent rows than before -> update the best linear independent rows

                    best_sizeR = count_R; % update the maximum number of linear independent rows found

                    best_R = R; % best set of linear independent rows
                    best_C = C; % best set of linear independent columns

                end


                if counter == maxIter && flag_Row == 0 % If it reaches the maximum number of iteration and we didn't find r linear independent rows

                    flag_Row = 2; % Update flag

                end

            end

            R = best_R; % Update Rows
            C = best_C; % Update Cols

            if flag_Row == 2 % If we didn't find r linear independent rows

                minSigmaFactor = abs(floor(log10(abs(min(svd(A(R,:),'econ')))))); % find min factor of singular values of A
                
                C = (1:n); % Use all columns

                % Start Greedy Light 

                sizeR = length(R);

                while sizeR < r

                    flag = 0;

                    sizeRows = length(R);
                    sizeRows = sizeRows+1;

                    for i = (1:m)

                        X = isempty(find(R==i));

                        if X==1
            
                            R_save = [R;i];

                            min_sigma = min(svd(A(R_save,C),'econ'));

                            if min_sigma > 10^(-minSigmaFactor) 

                                flag = 1;

                                R = [R;i];
                                break
                            end
                        end
                    end

                    if flag == 0
                       minSigmaFactor = minSigmaFactor + 1; 
                    end

                    sizeR = length(R);

                end % Finish Greedy Light
            end % end if we didn't find r linear independent rows
            
            % Should be selected the same LI columns as LI rows to be sym
            
            C = R; 
            
        else % if A is not symmetric
            
            counter = 0; % counter of random sets of cols

            flag_Row = 0; % to start while cycle

            C_tot = [1:n];  % All columns

            best_sizeR = -1; % Best size of linear independent rows

            while flag_Row == 0 % start while cycle

                counter = counter + 1; % update counter

                if counter > 1 % if counter > 1 -> reorder

                    rng('shuffle','twister'); % control random number generator

                    C_tot = (randperm(length(C_tot))); % reorder columns

                end

                % Penalized Identity matrix above A

                A_new = [zeros(r,n); A]; % new matrix A with the penalized identity submatrix above 
                
                penValue = penValue_Min;

                m1 = m + r; % Number of rows of the matrix A with the penalized identity submatrix above

                R  = [1:r]; % initial rows of the submatrix
                Rb = [r+1:m1]; % initial rows that aren't in the submatrix

                C = C_tot(1:r); % columns of the submatrix

                for i = (1:r)

                    id_R = R(i); % submatrix Row indice
                    id_C = C(i); % submatrix Col indice

                    A_new(id_R,id_C) = 1* 10^(-penValue); % replace 0 to penalized value to create an Linear Independent submatrix

                end

                countRi = 0; % counter of Linear Independent Rows founded

                flag_LSFI = 1; % to start the local search changing rows

                while flag_LSFI == 1 && countRi < r

                    flag_LSFI = 0;

                    Ar = A_new(R,C);

                    [L1,U1,P1] = lu(Ar');

                    i = 1;

                    while i <= m1-r-countRi

                        % LU Factoration   

                        b = P1 * A_new(Rb(i),C)';

                        y = L1\b;

                        alfa = U1\y;

                        % Changing Ar

                        local_alfa = find(abs(alfa)>1);

                        T = isempty(local_alfa);

                        if T == 0

                            [tamAlfa,~] = size(local_alfa);

                            flagAlfa = 0;

                            for iAlfa = (1:tamAlfa)

                                valueAlfa = local_alfa(iAlfa);

                                if R(valueAlfa) <= r %Alteracao v2, nao tinha R no valueAlfa

                                    local_alfa = valueAlfa;
                                    flagAlfa = 1;

                                    break
                                end
                            end

                            if flagAlfa == 0
                                local_alfa = local_alfa(1);

                            end

                            if R(local_alfa) <= r % Preference to change rows of Penalized Identity submatrix 

                                countRi  = countRi + 1;

                                R(local_alfa) = Rb(i);
                                Rb(i) = []; %Remove element from matrix

                                flag_LSFI = 1;

                                if countRi==r
                                    break
                                end
                            else

                                local_alfa = local_alfa(1);

                                el_save = R(local_alfa) ;

                                R(local_alfa) = Rb(i);

                                Rb(i) = el_save;

                            end

                            Ar = A_new(R,C);

                            [L1,U1,P1] = lu(Ar');
                        end
                        i = i + 1;
                    end % end while
                    
                    if flag_LSFI == 0 && countRi < r && penValue < penValue_Max

                        penValue = penValue + 1;

                        for i = (1:r)

                            id_R = R(i); % submatrix Row indice

                            if id_R <= r

                                id_C = C(i); % submatrix Col indice
                                A_new(id_R,id_C) = 1* 10^(-penValue); %update penalized value to create an Linear Independent submatrix

                            end
                        end
                        flag_LSFI = 1;
                    end
                    
                end % End of local search

                R_init = R; % copy R_init from R

                R = []; % Actual linear independent rows of A
                count_R = 0; % number of linear independent rows of A

                for j = (1:r)

                    Rj_temp = R_init(j) - r; % A_new has r more elements than A. If this indice is negative or 0, so it doesn't belongs to A

                    if Rj_temp > 0 % if it belongs to A

                        R = [R;Rj_temp]; % Add this row to the set of Linear Independent Rows

                        count_R = count_R + 1; % Update linear independent rows counter
                    end
                end

                if count_R == r % Can finish the algorithm, because it was found the r linear independent rows

                    flag_Row = 1; % Update flag

                end

                if count_R > best_sizeR % if we found more linear independent rows than before -> update the best linear independent rows

                    best_sizeR = count_R; % update the maximum number of linear independent rows found

                    best_R = R; % best set of linear independent rows
                    best_C = C; % best set of linear independent columns

                end


                if counter == maxIter && flag_Row == 0 % If it reaches the maximum number of iteration and we didn't find r linear independent rows

                    flag_Row = 2; % Update flag

                end

            end

            R = best_R; % Update Rows
            C = best_C; % Update Cols

            if flag_Row == 2 % If we didn't find r linear independent rows

                minSigmaFactor = abs(floor(log10(abs(min(svd(A(R,:),'econ')))))); % find min factor of singular values of A
                
                C = (1:n); % Use all columns

                % Start Greedy Light 

                sizeR = length(R);

                while sizeR < r

                    flag = 0;

                    sizeRows = length(R);
                    sizeRows = sizeRows+1;

                    for i = (1:m)

                        X = isempty(find(R==i));

                        if X==1
            
                            R_save = [R;i];

                            min_sigma = min(svd(A(R_save,C),'econ'));

                            if min_sigma > 10^(-minSigmaFactor) 

                                flag = 1;

                                R = [R;i];
                                break
                            end
            
                        end
                    end

                    if flag == 0
                       minSigmaFactor = minSigmaFactor + 1; 
                    end

                    sizeR = length(R);

                end % Finish Greedy Light

                % Penalized Identity matrix on the left of A

                A_new = [zeros(m,r) A];  % new matrix A with the penalized identity submatrix on the left
                
                penValue = penValue_Min;

                for i = (1:r)

                    id_R = R(i); % submatrix Row indice
                    id_C = C(i); % submatrix Col indice

                    A_new(id_R,id_C) = 1*10^(-penValue); % replace 0 to penalized value to create an Linear Independent submatrix

                end

                n1 = n + r; % number of columns of the matrix A with the penalized identity submatrix above

                C  = [1:r]; % initial columns of the submatrix

                Cb = [r+1:n1]; % initial rows that aren't in the submatrix

                countCi = 0; % counter of Linear Independent Rows founded

                flag_Col = 1; % flag for while cycle

                % Start local search to change columns

                while flag_Col == 1 && countCi < r

                    flag_Col = 0;

                    Ar = A_new(R,C);

                    % FOR COLUMNS

                    [L2,U2,P2] = lu(Ar);

                    i = 1;

                    while i <= n1-r-countCi

                        b = P2 * A_new( R,Cb(i));

                        y = L2\b;

                        alfa = U2\y;

                        % Changing Ar

                        local_alfa = find(abs(alfa)>1);

                        T = isempty(local_alfa);

                        if T == 0

                            [tamAlfa,~] = size(local_alfa);

                            flagAlfa = 0;

                            for iAlfa = (1:tamAlfa)

                                valueAlfa = local_alfa(iAlfa);

                                if C(valueAlfa) <= r %Alteracao v2, nao tinha C no valueAlfa

                                    local_alfa = valueAlfa;
                                    flagAlfa = 1;

                                    break
                                end
                            end

                            if flagAlfa == 0
                                local_alfa = local_alfa(1);

                            end

                            if C(local_alfa) <= r

                                countCi  = countCi + 1;

                                C(local_alfa) = Cb(i);
                                Cb(i) = []; %Remove element from matrix

                                flag_Col = 1;

                                if countCi==r
                                    break
                                end


                            else

                                local_alfa = local_alfa(1);

                                el_save = C(local_alfa) ;

                                C(local_alfa) = Cb(i);

                                Cb(i) = el_save;
                            end

                            Ar = A_new(R,C);

                            [L2,U2,P2] = lu(Ar);
                        end
                        i = i + 1;
                    end
                                    
                    if flag_Col == 0 && countCi < r && penValue < penValue_Max

                        penValue = penValue + 1;

                        for i = (1:r)

                            id_C = C(i); % submatrix Row indice

                            if id_C <= r

                                id_R = R(i); % submatrix Col indice
                                A_new(id_R,id_C) = 1* 10^(-penValue); %update penalized value to create an Linear Independent submatrix

                            end
                        end
                        flag_Col = 1;
                    end
                end % end local search

                C_init = C; % copy C_init from C

                C = []; % Actual linear independent columns of A
                count_C = 0; % number of linear independent rows of A

                for j = (1:r)

                    Cj_temp = C_init(j) - r; % A_new has r more elements than A. If this indice is negative or 0, so it doesn't belongs to A

                    if Cj_temp > 0 % if it belongs to A

                        C = [C;Cj_temp]; % Add this row to the set of Linear Independent Rows

                        count_C = count_C + 1; % Update linear independent rows counter
                    end
                end

                if count_C < r % if it wasn't found all r cols linear independent, flag_Col = 2
                    flag_Col = 2;
                end


                if flag_Col == 2 %Go to Greedy Light if we don't have r linear independent cols

                    minSigmaFactor = abs(floor(log10(abs(min(svd(A(R,C),'econ')))))); % find min factor of singular values of A
               
                    % Start Greedy Light 

                    sizeC = length(C);

                    while sizeC < r

                        flag = 0;

                        sizeColumns = length(C);
                        sizeColumns = sizeColumns+1;

                        for i = (1:n)

                            X = isempty(find(C==i));

                            if X==1
 
                                C_save = [C;i];

                                min_sigma = min(svd(A(R,C_save),'econ'));

                                if min_sigma > 10^(-minSigmaFactor)

                                    flag = 1;

                                    C = [C;i];

                                    break
                                end
            
                            end

                        end

                        if flag == 0
                           minSigmaFactor = minSigmaFactor + 1; 
                        end

                        sizeC = length(C);

                    end % end Greedy Light
                end % end if didn't find all r linear independent cols
            end % end if didn't find all r linear independent rows
        end % end if A is not symmetric  
    end % end if m> n or m < n or m = n
    end % enf if rank(A,10^-6) < r
    
    time = toc;
    
end
