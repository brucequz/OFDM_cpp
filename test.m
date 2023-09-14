
pilot_rx = [complex(0.24783945228802418,-0.030681970026679162)
  complex(0.78713337814658634, -1.0443158783729976)
  complex(0.72114804998840243,1.7728081627140173)
  complex(0.47485493354382857, -1.6047170785159688)
  complex(0.34060963642657555, 0.78879927773285785)
  complex(0.90881413679803424, -0.17666138679578569)
  complex(0.039373320654956218, 1.6244197719465174)
  complex(0.97945262541722122, -1.5189682516372347)
  complex(0.30452160546993967, 1.3554017713298585)
  complex(-0.75665266670899312, 0.69868124005465126)
  complex(0.39794924639135187, -0.89089164263055442)
  complex(-0.83429033839660382, 0.28822843081654448)
  complex(0.5287854019270688, -1.5574406215875121)
  complex(1.2641481194601822, -1.5698798739001871)
  complex(0.34266973005210638, 1.5963591612886154)
  complex(0.39893476181730736, -1.0901435937752773)]';

pilot_tx = [complex(1,1)
 complex(-1,1)
 complex(1,-1)
 complex(-1,-1)
 complex(1,1)
 complex(-1,1)
 complex(1,-1)
 complex(-1,-1)
 complex(1,1)
 complex(-1,1)
 complex(1,-1)
 complex(-1,-1)
 complex(1,1)
 complex(-1,1)
 complex(1,-1)
 complex(-1,-1)]';

H_est = mmse_channel_estimator_per_subcarrier(pilot_rx, pilot_tx, 0.5/16);
            

function H_estimated = mmse_channel_estimator_per_subcarrier(rx_pilot, tx_pilot, noise_var)
    num_subcarriers = size(tx_pilot, 2); % Assuming tx_pilot is a matrix with pilot symbols for each subcarrier
    H_estimated = zeros(1, num_subcarriers);

    for k = 1:num_subcarriers
        % Extract pilot symbols for the k-th subcarrier
        rx_pilot_k = rx_pilot(:, k);
        tx_pilot_k = tx_pilot(:, k);

        % Calculate the cross-correlation between received and transmitted pilots
        R_yx = tx_pilot_k' * rx_pilot_k;

        % Calculate the auto-correlation of transmitted pilots
        R_xx = tx_pilot_k' * tx_pilot_k;

        % Calculate the MMSE channel estimate for the k-th subcarrier
        H_estimated(k) = R_yx / (R_xx + noise_var);
    end
end

