%Ma trận kiểm tra
H = [1 1 0 1 0 0
     0 1 1 0 1 0
     1 0 0 0 1 1
     0 0 1 1 0 1];

%Từ mã truyền
c = [0 0 1 0 1 1];

% Từ mã nhận
r = [1 0 1 0 1 1];

y = r;

%Số lần lặp
maxiter = 20;
iter = 0;

%Đánh dấu giải mã chưa thành công
success = 0;

while iter < maxiter && ~success
    %Khởi tạo ma trận lỗi
    E = zeros(4, 6);
    for j = 1:4
        for i = 1:6
            if H(j, i) == 1
                %Xác định ma trận lỗi
                E(j, i) = mod(sum(y .* H(j, :)), 2);
            end
        end
    end

    % Tìm vị trí có nhiều lỗi nhất
    for i = 1:6
        M(i) = sum(E(:, i));
    end
    [~, index] = max(M);

    % Sửa lỗi
    if M(index) ~= 0
        y(index) = mod(y(index) + 1, 2);
    end

    % Kiểm tra sau khi đảo bit
    areErrorsPresent = check_errors(H, y);
    if areErrorsPresent == 0 % Không lỗi
        success = 1;
        disp("No error");
    else %Có lỗi
        disp("Still errors");
    end
    iter = iter + 1;
end

% Hàm kiểm tra lỗi
function res = check_errors(H, current_frame)
    syndrome = mod(H * current_frame', 2); 
    areErrors = any(syndrome); 
    res = areErrors;
end
