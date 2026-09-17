#include <iostream>
#include <fstream>
#include <vector>   // 新增 vector
#include <string>   // 新增 string
#include <iomanip>  // 為了 setw, setfill
#include <cmath>
// #include "SA.h"
#include "SA_v3.h"
using namespace std;

// 這些全域變數保持不變
bf_data_t out_data_array_real[5][NR][no_lines * no_lines_phi];
bf_data_t out_data_array_imag[5][NR][no_lines * no_lines_phi];
bf_data_t_64 out_data_array[NR][no_lines * no_lines_phi];
bf_data_t_64 SA_partial_image_array[5][NR][no_lines * no_lines_phi];
void output_topk(string filename, ap_uint<3> TX_idx){
    ofstream outfile(filename);
    if(!outfile.is_open()){
        cout << "file open error: " << filename << endl; // 稍微改良錯誤訊息以便除錯
        return;
    }
    bf_data_t top_real[K];
    bf_data_t top_imag[K];
    for(int i = 0; i < no_lines * no_lines_phi; i++){
        // sort out_data_array_real and out_data_array_imag to get top K
        for(int j = 0; j < NR; j++){
            for(int k = j; k < NR; k++){
                float power_j = (float)out_data_array_real[TX_idx][j][i] * (float)out_data_array_real[TX_idx][j][i] +
                                (float)out_data_array_imag[TX_idx][j][i] * (float)out_data_array_imag[TX_idx][j][i];
                float power_k = (float)out_data_array_real[TX_idx][k][i] * (float)out_data_array_real[TX_idx][k][i] +
                                (float)out_data_array_imag[TX_idx][k][i] * (float)out_data_array_imag[TX_idx][k][i];
                if(power_k > power_j){
                    // swap
                    bf_data_t temp_real = out_data_array_real[TX_idx][j][i];
                    bf_data_t temp_imag = out_data_array_imag[TX_idx][j][i];
                    out_data_array_real[TX_idx][j][i] = out_data_array_real[TX_idx][k][i];
                    out_data_array_imag[TX_idx][j][i] = out_data_array_imag[TX_idx][k][i];
                    out_data_array_real[TX_idx][k][i] = temp_real;
                    out_data_array_imag[TX_idx][k][i] = temp_imag;
                }
            }
        }
        // output top K
        for(int j = 0; j < K; j++){
            float real = (float)out_data_array_real[TX_idx][j][i];
            float imag = (float)out_data_array_imag[TX_idx][j][i];
            outfile << "real : " << real << ", imag: " << imag <<  " | ";
        }
        outfile << endl;
    }
    outfile.close();
}
void output_image(string filename, int TX_idx){
    ofstream outfile(filename);
    if(!outfile.is_open() ){
        cout << "file open error: " << filename << endl; // 稍微改良錯誤訊息以便除錯
        return;
    }
    for(int i = 0; i < no_lines * no_lines_phi; i++){
        for(int j = 0; j < NR; j++){
            float real = (float)out_data_array_real[TX_idx][j][i];
            float imag = (float)out_data_array_imag[TX_idx][j][i];
            outfile << (real * real + imag * imag) << " ";
        }
        outfile << endl;
    }
    
    outfile.close();
}

bool output_SA_image(string filename, string dataset_name){
    bool flag = 0;
    std::ostringstream debug_out, debug_out2;
    debug_out << "/home/hung52852/SA/testdata/" << dataset_name << "/SA_mismatch_debug.txt";
    debug_out2 << "/home/hung52852/SA/testdata/" << dataset_name << "/SA_add_record.txt";
    ofstream outfile(filename);
    ofstream debug_add_file(debug_out2.str());
    ofstream debug_file(debug_out.str());
    if(!outfile.is_open()){
        cout << "file open error: " << filename << endl;
        return true; 
    }
    for(int i = 0; i < no_lines * no_lines_phi; i++){
        for(int j = 0; j < NR; j++){
            // add the five frames together
            float sum_real = 0;
            float sum_imag = 0;
            int non_zero_count = 0;
            for(int k = 0; k < 5; k++){
                sum_real += (float)out_data_array_real[k][j][i];
                sum_imag += (float)out_data_array_imag[k][j][i];
                if(out_data_array_real[k][j][i] != 0 || out_data_array_imag[k][j][i] != 0)
                    non_zero_count += 1;
            }
            if(non_zero_count >= 1){
                debug_add_file << "depth: " << j << ", line: " << i / no_lines << " ,phi: " << i % no_lines << ": non-zero TX count = " << non_zero_count << endl;
                debug_add_file << "  Sum real = " << sum_real << ", Sum imag = " << sum_imag << endl;
                debug_add_file << "  Individual contributions: ";
                for(int k = 0; k < 5; k++){
                    debug_add_file << "(" << (float)out_data_array_real[k][j][i] << ", " << (float)out_data_array_imag[k][j][i] << ") ";
                }
                debug_add_file << endl;
            }
            float power = sum_real * sum_real + sum_imag * sum_imag;
            outfile << out_data_array[j][i] << " ";
            float diff = power - (float)out_data_array[j][i];
            float err = abs(diff / (float)out_data_array[j][i]);
            if(err > 0.05){ // 5% error threshold
                cout << "Mismatch at depth " << j << ", line :  " << i / no_lines << " , phi: " << i % no_lines << ": computed " << power << ", stored " << out_data_array[j][i] << endl;
                debug_file << "Mismatch at depth " << j << ", line " << i / no_lines << ", phi " << i % no_lines << ": computed " << power << ", stored " << out_data_array[j][i] << endl;
                flag = 1;
            }
        }
        outfile << endl;
    }
    outfile.close();
    debug_file.close();
    debug_add_file.close();
    return flag;
}

int main(){
    // 定義你要跑的所有資料夾名稱
    //vector<string> datasets = {"30_30_50","30_30_100", "30_30_150", "0_0_100", "0_0_150", "0_0_50"};
    //vector<string> datasets = {"0_0_100", "0_0_150", "0_0_50"};
    //vector<string> datasets = { "-30_-30_50", "-30_-30_100", "-30_-30_150", "-30_30_50", "-30_30_100", "-30_30_150" , "30_-30_50", "30_-30_100", "30_-30_150" };
    vector<string> datasets = { "0_0_100" };
    // 最外層迴圈：遍歷每個資料夾
    for(const string& dataset_name : datasets) {
        cout << "==================================================================================" << endl;
        cout << "         Processing Dataset: " << dataset_name << endl;
        cout << "==================================================================================" << endl;

        // 處理 5 個 TX
        for(int i = 0; i < 5; i++){
            std::ostringstream ss_real, ss_imag, out;
            bool reset;
            reset = (i == 0) ? true : false; // 只有第一個 TX 時重置 SA_image
            // 修改路徑：使用 dataset_name 變數
            ss_real << "/home/hung52852/SA/testdata/" << dataset_name << "/bb_tx"
                    << std::setw(2) << std::setfill('0') << (i+1)
                    << "_real.txt";

            ss_imag << "/home/hung52852/SA/testdata/" << dataset_name << "/bb_tx"
                    << std::setw(2) << std::setfill('0') << (i+1)
                    << "_imag.txt";

            ifstream infile_imag(ss_imag.str());
            ifstream infile_real(ss_real.str());

            if(!infile_imag.is_open() || !infile_real.is_open()){
                cout << "Input file open error for dataset " << dataset_name << " at TX " << i+1 << endl;
                continue; // 如果讀不到檔，跳過這個 TX 繼續下一個，而不是直接 return -1 結束程式
            }

            out << "/home/hung52852/SA/testdata/" << dataset_name << "/beamformed_output_tx"
                << std::setw(2) << std::setfill('0') << (i+1)
                << ".txt";
            
            ofstream outfile(out.str());

            bf_data_t_int rx_data_real[Data_dim * channel_num];
            bf_data_t_int rx_data_imag[Data_dim * channel_num];
            ap_uint<3> TX_idx = i;
            hls::stream<bf_data_t_64> out_data;
            //hls::stream<bf_data_t_64> empty_partial_image; // 空的 stream 傳入
            hls::stream<bf_data_t_64> SA_partial_image; // 用來接收部分影像的 stream
            for(int k = 0; k < Data_dim; k++){ // 注意：這裡原本變數名稱是 i，與外層迴圈衝突，建議改為 k 或其他
                for(int j = 0; j < channel_num; j++){
                    int idx = k * channel_num + j;
                    infile_real >> rx_data_real[idx];
                    infile_imag >> rx_data_imag[idx];
                }
            }
            // for(int i = 0; i < NR; i++){
            //     for(int j = 0; j < no_lines * no_lines_phi; j++){
            //         empty_partial_image.write(0); // 寫入空資料
            //     }
            // }
            top_model(reset, TX_idx, rx_data_real, rx_data_imag, out_data, SA_partial_image);

            // for(int k = 0; k < NR; k++){
            //     float sum = 0;
            //     for(int j = 0; j < no_lines * no_lines_phi; j++){
            //         bf_data_t_64 data_out = out_data.read();
            //         out_data_array[k][j] = data_out;
                    
            //         bf_data_t_64 partial_image_data = SA_partial_image.read();
            //         bf_data_t real;
            //         real.range(31,0) = partial_image_data.range(63,32);
            //         bf_data_t imag;
            //         imag.range(31,0) = partial_image_data.range(31,0);
            //         // 這裡存入全域陣列
            //         out_data_array_real[i][k][j] = real;
            //         out_data_array_imag[i][k][j] = imag;
                    
            //         sum += (float)(real * real + imag * imag);
            //     }
            //     outfile << sum << endl;
            // }
            
            float sum[NR] = {0}; // 將 sum 陣列移到外層迴圈外，並初始化為 0
            for (int j = 0; j < no_lines; j++) {
                for (int n = 0; n < no_lines_phi; n++) {
                    for (int k = 0; k < NR; k++) {
                        int idx = j * no_lines_phi + n;

                        bf_data_t_64 data_out;
                        out_data.read(data_out);
                        out_data_array[k][idx] = data_out;   // 這裡才對位

                        bf_data_t_64 partial_image_data;
                        partial_image_data = SA_partial_image.read();
                        bf_data_t real;
                        real.range(31,0) = partial_image_data.range(63,32);
                        bf_data_t imag;
                        imag.range(31,0) = partial_image_data.range(31,0);
                        // 這裡存入全域陣列
                        out_data_array_real[i][k][idx] = real;
                        out_data_array_imag[i][k][idx] = imag;  
                        sum[k] += (float)(real * real + imag * imag);
                    }
                }
            }
            std::cout << "out_data size = " << out_data.size() << "\n";
            std::cout << "SA_partial_image size = " << SA_partial_image.size() << "\n";
            for (int k = 0; k < NR; k++) {
                outfile << sum[k] << "\n";
            }
            outfile << std::endl;
            // cout << "Dataset: " << dataset_name << " | Test case " << i+1 << " completed." << endl;
            
            infile_imag.close();
            infile_real.close();
            outfile.close();
            std::cout << "=========================================\n";
            std::cout << "Completed TX " << (i+1) << " for dataset " << dataset_name << std::endl;
            std::cout << "Output written to: " << out.str() << std::endl;
            std::cout << "=========================================\n";
        }

        // Output individual images for this dataset
        for(int i = 0; i < 5; i++){
            std::ostringstream ss;
            ss << "/home/hung52852/SA/testdata/" << dataset_name << "/SA_image_tx"
               << std::setw(2) << std::setfill('0') << (i+1)
               << ".txt";
            output_image(ss.str(), i);
        }

        // Output the summed SA image for this dataset
        std::ostringstream ss_sum;
        ss_sum << "/home/hung52852/SA/testdata/" << dataset_name << "/SA_image_sum.txt";
        bool flag = output_SA_image(ss_sum.str(), dataset_name);
        if(flag == true){
            cout << "Mismatch found in summed SA image for dataset " << dataset_name << endl;
            //return flag;
        }
        for (int i = 0; i < 5; i++)
        {
            std::ostringstream ss_topk;
            ss_topk << "/home/hung52852/SA/testdata/" << dataset_name << "/SA_topk_tx"
                << std::setw(2) << std::setfill('0') << (i+1)
                << ".txt";
            output_topk(ss_topk.str(), i);
        }
        
        cout << "Done with " << dataset_name << endl << endl;
    }

    return 0;
}