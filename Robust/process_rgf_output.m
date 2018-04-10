for DATA = 3
    for NORM = [0 1]
        for EPS = [0.5 1]
            for UNC = [0.1 0.2]

                new_output_us = [];
                new_data_us = [];
                new_feedback_us = [];

                for INSTANCE = 1:50

                    filename = strcat('mat/2016-08-26_cluster_runs/output', num2str(DATA), ...
                    '_norm', num2str(NORM), '_eps', num2str(EPS), ...
                    '_unc', num2str(UNC), '_instance', num2str(INSTANCE), ...
                    '.mat');

                    if(exist(filename) ~= 2)
                        fprintf('File does not exist\n');
                    else
                        load(filename);
                        new_output_us{INSTANCE} = output_us{INSTANCE};
                        new_data_us{INSTANCE} = data_us{INSTANCE};
                        new_feedback_us{INSTANCE} = feedback_us{INSTANCE};
                    end

                end

                filename = strcat('mat/output', num2str(DATA), ...
                '_norm', num2str(NORM), '_eps', num2str(EPS), ...
                '_unc', num2str(UNC), '.mat');

                output_us = new_output_us;
                data_us = new_data_us;
                feedback_us = new_feedback_us;

                save(filename, 'output_us', 'data_us', 'feedback_us');

            end
        end 
    end
end
