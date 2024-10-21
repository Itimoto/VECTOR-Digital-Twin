%% Naive Time of Flight Estimation
% https://tns.thss.tsinghua.edu.cn/wst/docs/features/
function [tof_mat] = naive_tof(csi_data, bw)
    % naive_tof
    % Input:
    %   - csi_data is the CSI used for ranging; [T S A E]
    %       - ifft_point and bw are the IFFT and bandwidth parameters
    %   - bw is the channel bandwidth
    % Output:
    %   - tof_mat is the rough time-of-flight estimation result [T A]
    %
    % T indicates the number of CSI Packets
    % S indicates the number of subcarriers
    % A indicates the number of antennas (ie. the STS number in MIMO system)
    % L indicates the number of HT-LTFs in a single PPDU
    c = physconst('LightSpeed');
    [packet_num, subcarrier_num, antenna_num, extra_num] = size(csi_data);
    ifft_point = power(2, ceil(log2(subcarrier_num)));
    % Get CIR from each packet and each antenna by ifft(CFR);
    cir_sequence = zeros(packet_num, antenna_num, extra_num, ifft_point);

    for p = 1:packet_num
        for a = 1:antenna_num
            for e = 1:extra_num
                cir_sequence(p, a, e, :) = ifft(csi_data(p, :, a, e), ifft_point);
            end
        end
    end

    cir_sequence = squeeze(mean(cir_sequence, 4)); % [T ifft_point A]
    half_point = ifft_point / 2;
    half_sequence = cir_sequence(:, 1:half_point, :); % [T half_point A]
    peak_indices = zeros(packet_num, antenna_num); % [T A]

    for p = 1:packet_num
        for a = 1:antenna_num
            [~, peak_indices(p, a)] = max(half_sequence, [], 2);
        end
    end
    % Calculate ToF;
    tof_mat = peak_indices .* subcarrier_num ./ (ifft_point .* bw); % [T A]
end