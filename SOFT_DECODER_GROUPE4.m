% /!\ on admet code le code LDPC est régulier

% HARIRI, JERBI, LEVRARD, RUFFLE, SAMI 

% ON NOTE :
% c_node -> index j of c_node f_j
% v_node -> index i of v_node c_i
% r0(c_node, v_node) -> r_ji(0)
% q0(c_node, v_node) -> q_ij(0)
% q1 -> q_ij(1)
% Q0 -> Q_i(0)
% Q1 -> Q_i(1)

function [c_cor] = SOFT_DECODER_GROUPE0(c_ds_flip, H, P1_ds, MAX_ITER)
    c_cor = c_ds_flip; % au début le mot code corrigé est le mot code reçu
    r0 = zeros( size(H,1), size(H,2) ); % stocke les réponses r_ij(0) des c_node, les r_ij(1) ne sont pas stockés car r_ij(1) = 1 - r_ij(0) 
    q0 = zeros( size(H,1), size(H,2) ); % stocke les réponses q_ij(0) des v_node, les q_ij(1) ne sont pas stockés car q_ij(1) = 1 - q_ij(0) 

    % ETAPE 1
    % initialement les réponses q(0) sont les proba 1-P1_ds
    for c_node = 1:size(q0, 1)
        for v_node = 1:size(q0, 2)
            if H(c_node, v_node)
                q0(c_node, v_node) = 1 - P1_ds(v_node); 
            end
        end
    end

    for iter = 1:MAX_ITER
    %fprintf('ITER %i\n', iter);

        % ETAPE 2
        % calcul des réponses des c_node
        % (pour réaliser le produit de tous les q1 sauf pour le v_node considéré
        % on utilise 2 boucles for et on compare v_node et v_node_bis
        % si ils sont différents on continu le produit)
        for c_node = 1:size(H,1) % for j (j est un c_node)
            for v_node = 1:size(H,2) % for i (i est un v_node)
                if H(c_node, v_node) % les deux nodes sont liés
                    r0(c_node, v_node) = 1;
                    for v_node_bis = 1:size(H,2)
                        if v_node_bis ~= v_node && H(c_node, v_node_bis)
                            r0(c_node, v_node) = r0(c_node, v_node)*(1 - 2*( 1-q0(c_node, v_node_bis) )); % on utilise q1=1-q0
                            %fprintf('c_node=%i  | v_node=%i | v_node_bis=%i | q0=%d | r0=%d\n', c_node, v_node, v_node_bis, q0(v_node_bis), r0(c_node, v_node));
                        end
                    end
                    r0(c_node, v_node) = 0.5 + 0.5*r0(c_node, v_node); % on ajoute les composantes indépendantes de i et j
                    %fprintf('r0=%i\n', r0(c_node, v_node));                    
                end
            end
        end

        % ETAPE 3
        % - calcul des réponses des v_node
        % - calcul de Q0 et Q1
        % - choix du bit
        % - check de la parité pour terminer le programme
        for v_node = 1:size(H,2)        
            Q0 = 1; % Q0 et Q1 sont propre à un v_node
            Q1 = 1; % on les reset à chaque v_node car pas besoin de les stocker ; ils ne servent qu'à choisir le bit le plus probable pour un v_node
            for c_node = 1:size(H,1)
                if H(c_node, v_node)
                    q0(c_node, v_node) = 1;
                    q1 = 1; % variable temporaire car ne sert qu'à trouver K pour esnuite calculer q0
                    Q0 = Q0*r0(c_node, v_node);
                    Q1 = Q1*(1-r0(c_node, v_node));
                    for c_node_bis = 1:size(H,1) % même méthode qu'à l'ETAPE 2 : 2 boucles for pour le produit
                        if c_node_bis ~= c_node && H(c_node_bis, v_node)
                            q0(c_node, v_node) = q0(c_node, v_node)*r0(c_node_bis, v_node);
                            q1 = q1*(1-r0(c_node_bis, v_node));
                            %fprintf('v_node=%i | c_node=%i | c_node_bis=%i | r0=%d | q0=%d\n', v_node, c_node, c_node_bis, r0(c_node_bis, v_node),q0(v_node))
                        end
                    end
                    % multiplication par la proba (voir formule)
                    q1 = P1_ds(v_node)*q1;
                    q0(c_node, v_node) = (1-P1_ds(v_node))*q0(c_node, v_node);

                    K = 1/(q0(c_node, v_node)+q1); % calcul de K

                    q0(c_node, v_node) = K*q0(c_node, v_node);
                end
            end
            
            % calculs Q0 et Q1
            Q1 = P1_ds(v_node)*Q1;
            Q0 = (1-P1_ds(v_node))*Q0;

            K = 1/(Q0 + Q1); % calcul de K_ij

            Q0 = K*Q0;
            Q1 = K*Q1;
            
            % choix du bit le plus probable
            if Q1 > Q0
                c_cor(v_node) = 1;
            else
                c_cor(v_node) = 0;
            end
        end
        
        % vérification de la parité pour chaque c_node
        % si vérifiée, on peut stoper le programme
        all_check = true;
        for c_node = 1:size(H,1)
            sum = 0;
            for v_node = 1:size(H,2)
                if H(c_node, v_node)
                    sum = sum + c_cor(v_node); % pour vérifier la parité on somme les bits portés la les v_nodes lié au c_node
                end
            end
            if mod(sum, 2) % => true si vaut 1 => nombre impaire de 1 => pas de parité
                all_check = false;
                break;
            end
        end
    
        if all_check
            %disp('ALL CHECK');
            break; % on stope le programme
        end
    end
end