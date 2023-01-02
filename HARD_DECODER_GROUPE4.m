% /!\ on admet code le code LDPC est régulier

% HARIRI, JERBI, LEVRARD, RUFFLE, SAMI 

function [c_cor] = HARD_DECODER_GROUPE0(c, H, MAX_ITER)

c_cor = c; % au début le mot code corrigé est le mot code reçu

% calcul de w_c
w_c = 0; % nombre de 1 par colonne dans H => nb de c_node lié à chaque v_node
for row = 1:size(H,1)
    if H(row, 1)
        w_c = w_c + 1;
    end
end


% début des itérations de décodage
for iter = 1:MAX_ITER

    %fprintf('ITER %i\n', iter);
    all_check = true; % flag "parité satisafaite pour chaque c_node"

    % CALCUL DE LA PARITÉ ET DE LA CORRECTION (ÉTAPE 2)
    % (en réalité on calcul la correction pour chaque v_node lié au c_node
    % et on la compare au bit porté le v_node, si c'est le même bit on en
    % déduit que la parité est vérifiée)

    c_node_sending = zeros( size(H,1) , size(H,2) ); % matrice qui socke les réponses des c_nodes, même dim que H pour simplifier l'accès aux cellules
    for c_node = 1:size(H,1) % size(H,1) est le nombre de ligne de H donc le nombre de c_node
        %   fprintf('c_node %i\n', c_node)
       
        % pour vérifier la parité du c_node on somme les bits portés par les v_node qui lui sont liés
        for v_node = 1:size(H,2) % on considère un v_node
            sum = 0;
            if H(c_node, v_node) % si ce v_node est lié au c_node
                for v_node_bis = 1:size(H,2) % pour calculer la correction on doit sommer les bit des autres v_node
                    if v_node_bis ~= v_node & H(c_node, v_node_bis) % si c'est un autre v_node que celui qu'on considère et si il est bien lié au c_node également
                        sum = sum + c_cor(v_node_bis); % on ajoute la valeur du v_node
                    end
                end
                
                cor = mod(sum, 2); % calcul du bit qui permet de vérifier la parité avec les autres bits
                c_node_sending(c_node, v_node) = cor; % stockage dans la matrice
                if cor ~= c_cor(v_node) % si différent => la parité n'etait pas vérifiée avant qu'on corrige
                    all_check = false;
                    
                end
            end
        end
        %   disp(c_node_sending);
    end

    % CHOIX DU BIT LE PLUS PROBABLE (ÉTAPE 3)
    % on réalise un majority vote parmi les différentes sources : les
    % réponses des c_node lié à un v_node + le bit porté la le v_node
    for v_node = 1:size(c_node_sending, 2)

        %fprintf('=> V_NODE %i\n', v_node);
        sum = c_cor(v_node); % somme initialisée avec valeur du bit du v_node
        %fprintf('init sum = %i\n', sum);

        for c_node = 1:size(c_node_sending, 1)
            if H(c_node, v_node) % si le c_node est lié au v_node
                sum = sum + c_node_sending(c_node,v_node); % on ajoute sa valeur
            end
        end
        %fprintf('sum = %i\n', sum);

        if sum >= (w_c+1)/2 % w_c+1 est le nombre de source
            c_cor(v_node) = 1;
            %disp('choose 1');
        else
            c_cor(v_node) = 0;
            %disp('choose 0');
        end
    end

    if all_check
        %disp('ALL CHECK')
        break;
    end
end
end