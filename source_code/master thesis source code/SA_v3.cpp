
#include <iostream>
#include <fstream>
#include "SA_v3.h"
#include "hls_math.h"
#include "sin_cos_value.h"
using namespace std;
/*
2026 / 1 /6 
for sparse memory access pattern, we need to change the loop order to
for j = 0 to no_lines_phi
    for k = 0 to no_lines
        for i = 0 to NR
This is to ensure that for every scan line we can find the top 30 values among all the scan points in that line.
The update of SA image is to see if the new and old data also belong to the SA image.
If the SA contain the new data, we add it. If the SA contain the old data, we minus it.
If both new and old data belong to SA image, we do add and minus accordingly.
If neither new nor old data belong to SA image, we do nothing.
*/ 
////////////////////////////////
//         top function       //
////////////////////////////////
//static ap_int<24> SA_image[no_lines][no_lines_phi][60] = {0}; // final SA image 
//static ap_int<9>  SA_image_coord[no_lines][no_lines_phi][60] = {0}; // final SA image coordinate
//static ap_uint<5> SA_mem_ptr[no_lines_phi][no_lines] = {0};
static sparse_mem_t partial_image[5][no_lines][no_lines_phi][K] ; // 5 TX, 30 partial sums per TX
static ap_uint<9>  partial_image_coord[5][no_lines][no_lines_phi][K] ; // 5 TX, 30 partial sums per TX
static ap_uint<284> partial_image_present_mask[5][no_lines][no_lines_phi]; // this is used to store the coord of the K partial sums in one ap_uint, it can replace the partial_image_coord 4D array to save memory
static ap_uint<54>  partial_image_valid_bit[5][no_lines][no_lines_phi]; // this is used to indicate whether the corresponding partial sum is valid or not
static ap_uint<9>   partial_image_min_coord[5][no_lines][no_lines_phi]; // this is used to store the min coord in the K partial sums
//static ap_uint<5> partial_mem_ptr[no_lines_phi][no_lines] = {0};
void top_model(bool reset, ap_uint<3> TX_idx, bf_data_t_int rx_data_real[Data_dim * channel_num], bf_data_t_int rx_data_imag[Data_dim * channel_num], 
                hls::stream<bf_data_t_64> &out_data, hls::stream<bf_data_t_64> &SA_partial_image){
    #pragma HLS INTERFACE axis port = out_data
    #pragma HLS INTERFACE s_axilite port = TX_idx bundle = control
    #pragma HLS INTERFACE m_axi port = rx_data_real bundle = raw_data_real
    #pragma HLS INTERFACE m_axi port = rx_data_imag bundle = raw_data_imag
    #pragma HLS INTERFACE s_axilite port = return bundle = control
    #pragma HLS array_partition variable = partial_image dim = 1 complete
    #pragma HLS array_partition variable = partial_image_present_mask dim = 1 complete
    // #pragma HLS array_partition variable = partial_image_coord dim = 1 complete
    // #pragma HLS array_partition variable = partial_image dim = 4 complete
    // #pragma HLS array_partition variable = partial_image_coord dim = 4 complete
    bf_data_t_int my_real_cache[Data_dim][channel_num];
    bf_data_t_int my_imag_cache[Data_dim][channel_num];
    // if(reset){
    //     #ifndef __SYNTHESIS__
    //     cout << "Reset SA_image_int to 0" << endl;
    //     #endif
    //     for(int TX_idx = 0; TX_idx < 5; TX_idx++){
    //         for(int i = 0; i < no_lines; i++){
    //             for(int j = 0; j < no_lines_phi; j++){
    //                 for(int k = 0; k < K; k++){
    //                     partial_image[TX_idx][i][j][k] = 0;
    //                     partial_image_coord[TX_idx][i][j][k] = 0;
    //                 }
    //             }
    //         }
    //     }
    // }
    load_local_cache(rx_data_real, my_real_cache);
    load_local_cache(rx_data_imag, my_imag_cache);
    // #ifndef __SYNTHESIS__
    //     for(int i = 0; i < channel_num; i++){
    //         for(int j = 0; j < Data_dim; j++){
    //             if (j > 1950 && j < Data_dim)
    //                 cout << "my_real_cache[" << j << "][" << i << "] = " << my_real_cache[j][i] << endl;
    //         }
    //     }
    // #endif
    SA_dataflow_top(my_real_cache, my_imag_cache, TX_idx, out_data, SA_partial_image);
    
}
void load_local_cache(bf_data_t_int rx_data[Data_dim * channel_num], bf_data_t_int local_cache_mem[Data_dim][channel_num]){
    #pragma HLS inline off
    for(int i = 0; i < Data_dim; i++){
        #pragma HLS pipeline II = 1
        for(int j = 0; j < channel_num; j++){
            local_cache_mem[i][j] = rx_data[i * channel_num + j];
            // #ifndef __SYNTHESIS__
            // if (i > 1450 && i < 1520)
            //     cout << "local_cache_mem[" << i << "][" << j << "] = " << local_cache_mem[i][j] << endl;
            // #endif
        }
    }
}
// void SA_dataflow_top(bf_data_t_int rx_data_real[Data_dim][channel_num], bf_data_t_int rx_data_imag[Data_dim][channel_num], ap_uint<3> TX_idx,
//                      hls::stream<bf_data_t_64> &out_data){
void SA_dataflow_top(bf_data_t_int rx_data_real[Data_dim][channel_num], bf_data_t_int rx_data_imag[Data_dim][channel_num], ap_uint<3> TX_idx, 
                     hls::stream<bf_data_t_64> &out_data, hls::stream<bf_data_t_64> &SA_partial_image){
    #pragma HLS dataflow
    hls::stream<trig_t> cos_phi;
    hls::stream<trig_t> sin_phi;
    hls::stream<trig_t> cos_th;
    hls::stream<trig_t> sin_th;
    hls::stream<inter_t> T_rd[channel_num];
    hls::stream<range_t> range;
    hls::stream<cos_sin_t> cos_stream[channel_num];
    hls::stream<cos_sin_t> sin_stream[channel_num];
    hls::stream<index_t> index_stream[channel_num];
    hls::stream<bf_data_t_64> partial_beamformed_stream;
    hls::stream<bf_data_t_64> partial_beamformed_stream2;
    #pragma HLS array_partition variable = cos_stream type = complete
    #pragma HLS array_partition variable = sin_stream type = complete
    #pragma HLS array_partition variable = index_stream type = complete
    #pragma HLS stream variable = partial_beamformed_stream 
    
    SA_beamform_seq_gen(cos_phi, sin_phi, cos_th, sin_th, T_rd, range);
    delay_approx(cos_phi, sin_phi, cos_th, sin_th, T_rd, range, TX_idx, index_stream, cos_stream, sin_stream);
    beamformer_CF(cos_stream, sin_stream, index_stream, rx_data_real, rx_data_imag, partial_beamformed_stream);
    //partial_beamformed_stream_U memory usage is crazy high, need to optimize
    // SA_sparse_mem_input_topK(partial_beamformed_stream, TX_idx, partial_beamformed_stream2);
    vector_adder(partial_beamformed_stream, TX_idx, out_data, SA_partial_image);
    #ifndef __SYNTHESIS__
        std::cout << "partial_beamformed_stream size = " << partial_beamformed_stream.size() << std::endl;
        std::cout << "partial_beamformed_stream2 size = " << partial_beamformed_stream2.size() << std::endl;
    #endif
}
/////////////////////////////////////
// Beamforming sequence generator  //
/////////////////////////////////////

void SA_beamform_seq_gen(hls::stream<trig_t> &cos_phi, hls::stream<trig_t> &sin_phi, hls::stream<trig_t> &cos_th,
                        hls::stream<trig_t> &sin_th, hls::stream<inter_t> T_rd[channel_num], hls::stream<range_t> &range) {
    
    seq_loop2:for(int j = 0; j < no_lines_phi; j++){
        seq_loop3: for(int k = 0; k < no_lines; k++){
            seq_loop: for(int i = 0; i < NR; i++){
                #pragma HLS PIPELINE II = 1
                //if((j == 34 && k == 34) || (j == 54 && k == 54)){ // we want 0 degree and 30 degree
                    cos_phi.write(cos_phi_lut[j]);
                    sin_phi.write(sin_phi_lut[j]);
                    cos_th.write(cos_th_lut[k]);
                    sin_th.write(sin_th_lut[k]);
                    range.write(range_lut[i]);
                    seq_ch_loop:for(int l = 0; l < channel_num; l++){
                        #pragma HLS UNROLL
                        T_rd[l].write(T_rd_lut[l][i]);
                    }
                //}
                // store or process the read values as needed
            }
        }
    }
}
void delay_approx(
/*input*/    hls::stream<trig_t> &cos_phi, hls::stream<trig_t> &sin_phi, hls::stream<trig_t> &cos_th, hls::stream<trig_t> &sin_th, hls::stream<inter_t> T_rd[channel_num], hls::stream<range_t> &range, ap_uint<3> TX_idx,
/*output*/   hls::stream<index_t> index[channel_num], hls::stream<cos_sin_t> cos_stream[channel_num], hls::stream<cos_sin_t> sin_stream[channel_num]){
    
    seq_loop2:for(int j = 0; j < no_lines_phi; j++){
        seq_loop3: for(int k = 0; k < no_lines; k++){
            seq_loop: for(int i = 0; i < NR; i++){
                #pragma HLS PIPELINE II = 1
                //if((j == 34 && k == 34) || (j == 54 && k == 54)){ // we want 0 degree and 30 degree
                    // TX_dist 
                    inter_t x_delta, y_delta, z_delta;
                    range_t range_val;
                    trig_t cos_phi_val, sin_phi_val, cos_th_val, sin_th_val;
                    range_val = range.read();
                    cos_phi_val = cos_phi.read();
                    sin_phi_val = sin_phi.read();
                    cos_th_val = cos_th.read();
                    sin_th_val = sin_th.read();
                    // Calculate deltas
                
                    x_delta = range_val * cos_phi_val * sin_th_val - TX_x_lut[TX_idx];
                    y_delta = range_val * sin_phi_val  - TX_y_lut[TX_idx];
                    z_delta = range_val * cos_phi_val * cos_th_val;

                    inter_t_sq TX_dist_sq = x_delta * x_delta + y_delta * y_delta + z_delta * z_delta;
                    inter_t TX_dist;
                    if(TX_idx == 0)
                        TX_dist = range_val;
                    else
                        TX_dist = hls::sqrt(TX_dist_sq);
                    // TX_dist = hls::sqrt(TX_dist_sq);
                    // if(i > NR - 10){
                    //     std::cout << "x_delta: " << x_delta << ", y_delta: " << y_delta << ", z_delta: " << z_delta << std::endl;
                    //     std::cout << "TX_dist: " << TX_dist << std::endl;
                    // }
                    // RX_dist and total_dist
                    channel_loop:
                    for(int l = 0; l < 64; l++){
                        #pragma HLS UNROLL
                        inter_t T_rd_temp = T_rd[l].read();
                        inter_t one_way_dist = T_rd_temp - chx[l] * cos_phi_val * sin_th_val - chy[l] * sin_phi_val;
                        accum_t total_dist = TX_dist + one_way_dist;
                        // Calculate delay
                        total_delay_t delay = total_dist * inv_c;
                        index_t index_temp = delay * fs;
                        // if (index_temp > (Data_dim - 1))
                        //     valid_stream[l].write(false);
                        // else 
                        //     valid_stream[l].write(true);
                        if(index_temp > (Data_dim - 1))
                            index_temp = Data_dim - 1;
                        // 2. 計算相位 phi = delay * fs * 2pi
                        
                        // phase_t phase = delay * 80000 * PI; // in radius
                        phase_t phase = delay * 40000; // normalize to angle
                        // 3. 截取小數部分 (0~1 turn) 送入 LUT
                        // 強制轉型為 <16, 0>，HLS 會自動捨棄整數部分 (Modulo 1.0)
                        // 這樣就完成了歸一化，且不需要任何額外邏輯
                        norm_phase_t theta_norm = phase;
                        cos_sin_t cos_val = approx_cos(theta_norm);
                        cos_sin_t sin_val = approx_sin(theta_norm);
                        cos_stream[l].write(cos_val);
                        sin_stream[l].write(sin_val);
                        index[l].write(index_temp);
                    }
                //}
                // store or process the read values as needed
            }
        }
        #ifndef __SYNTHESIS__
            std::cout << "Delay approximation for frame " << j << " completed." << std::endl;
        #endif
    }
}
trig_t approx_sin(norm_phase_t theta) {
    //#pragma HLS PIPELINE II=1
    #pragma HLS INLINE

    // 1. 取出象限資訊 (最高 2 bits)
    // theta 為 16 bits，bit 15-14 為象限
    // 00: Q1, 01: Q2, 10: Q3, 11: Q4
    ap_uint<2> quadrant = theta.range(15, 14);
    
    // 2. 取出查表 Index (接下來的 10 bits)
    // LUT_SIZE = 1024 = 2^10
    // bit 13-4 為 index
    ap_uint<10> index_raw = theta.range(13, 4);
    
    ap_uint<4>  index_frac = theta.range(3, 0); // 小數部分 (4 bits)
    // 3. 處理對稱性 (Symmetry)
    // Q1 (00) & Q3 (10): 波形上坡 -> 直接查表
    // Q2 (01) & Q4 (11): 波形下坡 -> 反向查表 (1023 - index)
    // 利用 quadrant[0] (LSB) 來判斷是否需要反轉
    // 在二進位中，(Max - Index) 等同於 bitwise NOT (~)
    ap_uint<10> index_eff, index_next;
    if (quadrant[0] == 1) { 
        //index_eff = ~index_raw;
        index_eff = (LUT_size - 1) - index_raw;
    } else {
        //index_eff = index_raw;
        index_eff = index_raw;
    }
    if(index_eff == (LUT_size - 1)){
        index_next = index_eff; // 避免溢位
    } else {
        index_next = index_eff + 1;
    }
    // 4. 查表
    // 您的 LUT 存的是 double/float，讀出來後轉為 fixed point
    trig_t val_0 = (trig_t) sin_cos_lut[index_eff];
    trig_t val_1 = (trig_t) sin_cos_lut[index_next];
    ap_fixed<32, 4> acc = val_0 + (val_1 - val_0) * ((ap_fixed<32, 4>)index_frac / 16); // 4 bits fraction -> divide by 16
    trig_t val = (trig_t) acc;  
    // 5. 處理正負號 (Sign)
    // Q1 (00) & Q2 (01): 正 (Upper half)
    // Q3 (10) & Q4 (11): 負 (Lower half)
    // 利用 quadrant[1] (MSB) 來判斷正負
    trig_t result;
    if (quadrant[1] == 1) { 
        result = -val;
    } else {
        result = val;
    }

    return result;
}
trig_t approx_cos(norm_phase_t theta) {
    //#pragma HLS PIPELINE II=1
    #pragma HLS INLINE
    // Cos(x) = Sin(x + pi/2)
    // 在歸一化系統 (0~1) 中，pi/2 = 0.25
    // 0.25 在 ap_ufixed<16, 0> 中等於 0b01000000... (次高位為1)
    // 我們可以利用無號數溢位特性直接相加
    
    // 方法 1: 直接加常數 (HLS 會優化成位元操作)
    // ap_ufixed<16, 0> offset = 0.25;
    // norm_phase_t theta_shifted = theta + offset;
    
    // 方法 2: 位元操作 (Toggle MSB-1)
    norm_phase_t theta_shifted = theta;
    theta_shifted[14] = ~theta[14]; // Toggle bit 14 (相當於 +0.25 or -0.75)

    return approx_sin(theta_shifted);
}
/////////////////////////////////////
//          Beamformer             //
/////////////////////////////////////
void beamformer_CF(hls::stream<cos_sin_t> cos_stream[channel_num],
                    hls::stream<cos_sin_t> sin_stream[channel_num],
                    hls::stream<index_t> idx_stream[channel_num],
                    bf_data_t_int rx_data_real[Data_dim][channel_num],
                    bf_data_t_int rx_data_imag[Data_dim][channel_num],
                    hls::stream<bf_data_t_64> &out_data) {
    #pragma HLS array_partition variable = rx_data_real dim = 2 type = complete
    #pragma HLS array_partition variable = rx_data_imag dim = 2 type = complete
    #pragma HLS array_partition variable = cos_stream type = complete
    #pragma HLS array_partition variable = sin_stream type = complete
    #pragma HLS array_partition variable = idx_stream type = complete
    bf_data_t_demod sample_real[channel_num];
    bf_data_t_demod sample_imag[channel_num];
    bf_data_t_demod sample_real_abs[channel_num];
    bf_data_t_demod sample_imag_abs[channel_num];
    #pragma HLS array_partition variable = sample_real type = complete
    #pragma HLS array_partition variable = sample_imag type = complete
    #ifndef __SYNTHESIS__
        ofstream demod_real("/home/hung52852/SA/demod_real.txt");
        ofstream demod_imag("/home/hung52852/SA/demod_imag.txt");
        std::cout << "Start beamformer_CF processing..." << std::endl;
        int counter = 0;
        index_t debug_index;
    #endif
    //while(true){
    
    beamformer_loop2:
    for(int i = 0; i < no_lines; i++){
        beamformer_loop3:
        for(int j = 0; j < no_lines_phi; j++){
            beamformer_loop:
            for(int k = 0; k < NR; k++){
            #pragma HLS pipeline II = 1
                for (int ch = 0; ch < channel_num; ch++) {
                    #pragma HLS unroll
                    cos_sin_t cos_phase = cos_stream[ch].read();
                    cos_sin_t sin_phase = sin_stream[ch].read();
                    index_t index_delay = idx_stream[ch].read();
                    
                    sample_real[ch] = rx_data_real[index_delay][ch] * cos_phase - rx_data_imag[index_delay][ch] * sin_phase;
                    sample_imag[ch] = rx_data_real[index_delay][ch] * sin_phase + rx_data_imag[index_delay][ch] * cos_phase;
                    sample_real_abs[ch] = (sample_real[ch] < 0) ? (bf_data_t_demod)-sample_real[ch] : sample_real[ch]; // mind the overflow here -128 to 127
                    sample_imag_abs[ch] = (sample_imag[ch] < 0) ? (bf_data_t_demod)-sample_imag[ch] : sample_imag[ch];
                    
                    
                    #ifndef __SYNTHESIS__
                        // demod_real << "Channel " << ch << ", Index " << index_delay << ": " << sample_real[ch] << std::endl;
                        // demod_imag << "Channel " << ch << ", Index " << index_delay << ": " << sample_imag[ch] << std::endl;
                        // if (counter / (no_lines * no_lines_phi) > (NR - 10)){
                        //     demod_real << "Channel " << ch << ", Index " << index_delay << ": " << sample_real[ch] << std::endl;
                        //     demod_imag << "Channel " << ch << ", Index " << index_delay << ": " << sample_imag[ch] << std::endl;
                        //     std::cout << "Channel " << ch << ", Index " << index_delay << ": " << sample_real[ch] << std::endl;
                        //     std::cout << "Channel " << ch << ", Index " << index_delay << ": " << sample_imag[ch] << std::endl;
                        // }
                    #endif
                }
                #ifndef __SYNTHESIS__
                    if (counter / (no_lines * no_lines_phi) > (NR - 10) && counter % (no_lines * no_lines_phi) == 0){
                        demod_real << counter / (no_lines * no_lines_phi) << " Channel " << 32 << ", Index " << debug_index << ": " << sample_real[32] << std::endl;
                        demod_imag << counter / (no_lines * no_lines_phi) << " Channel " << 32 << ", Index " << debug_index << ": " << sample_imag[32] << std::endl;
                    }
                #endif
                bf_data_t real_part = sum_pipeline(sample_real);
                bf_data_t imag_part = sum_pipeline(sample_imag);
                //bf_data_t CF_denom_real = sum_pipeline(sample_real_abs);
                //bf_data_t CF_denom_imag = sum_pipeline(sample_imag_abs);
                bf_data_t CF_denom_real = 0;
                bf_data_t CF_denom_imag = 0;
                for(int ch = 0; ch < channel_num; ch++){
                    #pragma HLS unroll
                    CF_denom_real += sample_real[ch] * sample_real[ch];
                    CF_denom_imag += sample_imag[ch] * sample_imag[ch];
                }
                bf_data_t CF_denom = CF_denom_real + CF_denom_imag ;
                
                //bf_data_t real_part_abs = (real_part < 0) ? (bf_data_t)-real_part : real_part;
                //bf_data_t imag_part_abs = (imag_part < 0) ? (bf_data_t)-imag_part : imag_part;
                bf_data_t CF;
                if(CF_denom != 0){
                    CF = (real_part * real_part + imag_part * imag_part) / CF_denom;
                }
                else{
                    CF = 0;
                }
                // the summation of SA should be in complex form
                // output the concade of real part and imag part
                //bf_data_t_64 output = (real_part * real_part + imag_part * imag_part) * CF * CF;
                bf_data_t_64 output;
                bf_data_t real_final, imag_final;
                real_final = real_part * CF; // surpress sidelobe , mind if it let the mainlobe also drop too much
                imag_final = imag_part * CF;
                //real_final = real_part; // test if CF affect the beampattern
                //imag_final = imag_part;
                output.range(63,32) = real_final.range();
                output.range(31,0)  = imag_final.range();
                out_data.write(output); 
                // bool end_signal = end_stream.read();
                #ifndef __SYNTHESIS__
                    counter ++;    
                    if (counter % (no_lines * no_lines_phi) == 0){
                        std::cout << "start processing " << counter / (no_lines * no_lines_phi) << " frames." << std::endl;
                    }
                #endif
                // if (end_signal == true) {
                //     break;
                //     #ifndef __SYNTHESIS__
                //         std::cout << "Beamformer_CF processing completed." << std::endl;
                //     #endif
                // }
        
    // }
            }
        }
    }
    #ifndef __SYNTHESIS__
        demod_real.close();
        demod_imag.close();
    #endif
}
bf_data_t sum_pipeline(bf_data_t_demod to_sum_element[channel_num]){
#pragma HLS inline off
    bf_data_t sum[channel_num] = {0};
    sum[0] = to_sum_element[0];
	#pragma HLS array_partition variable=sum dim = 0 type = complete
	#pragma HLS ARRAY_PARTITION variable=to_sum_element dim = 0 type=complete
    #pragma HLS pipeline II = 1
    //ofstream debug_sum("/home/hung52852/SA/debug_sum.txt");
    for(int i = 1; i < channel_num; i++){
		#pragma HLS unroll
        sum[i] = sum[i - 1] + to_sum_element[i];
        //debug_sum << "Intermediate sum at " << i << ": " << sum[i] << std::endl;
    }
    return sum[channel_num - 1];
}
///////////////////////////////////////////
// vector adder
///////////////////////////////////////////
// II = 2 need to be improved
void vector_adder(hls::stream<bf_data_t_64> &TX_partial_new_image, ap_uint<3> TX_idx, hls::stream<bf_data_t_64> &out_data, hls::stream<bf_data_t_64> &SA_partial_image){
    //ofstream debug_vector_adder("/home/hung52852/SA/debug_vector_adder.txt");
    #ifndef __SYNTHESIS__
        std::cout << "Start vector_adder processing..." << std::endl;
        std::ofstream debug_vector_adder("/home/hung52852/SA/debug_vector_adder_TX"
                                    + std::to_string(TX_idx)
                                    + ".txt");
        std::ostringstream sparse_mem;
        sparse_mem << "/home/hung52852/SA/testdata/debug_sparse_mem_TX"
                    << std::setw(2) << std::setfill('0') << (TX_idx+1)
                    << ".txt";
        std::ofstream debug_topk(sparse_mem.str());
    #endif
    vector_adder_loop1:for(int j = 0; j < no_lines; j++){
        vector_adder_loop2:for(int k = 0; k < no_lines_phi; k++){
            ap_uint<9>  top_coord[K] = {0};
            ap_uint<54> present_mask = 0;
            ap_uint<9>  min_coord = 0;
            ap_uint<6> valid_cnt = 0;
            sparse_mem_complex_t top_real[K] = {0};
            sparse_mem_complex_t top_imag[K] = {0};
            ap_uint<6> rptr[4] = {0}; // read pointer for the 4 other TX
            vector_adder_loop3:for(int i = 0; i < NR; i++){
                #pragma HLS PIPELINE II = 1
                bf_data_t_64 new_val = TX_partial_new_image.read();
                
                bf_data_t new_real, new_imag;
                new_real.range() = new_val.range(63,32);
                new_imag.range() = new_val.range(31,0);
                // ap_fixed<20,16> new_real_12;
                sparse_mem_complex_t new_real_12;
                new_real_12 = new_real.range(31,20);
                //new_real_12 = static_cast<ap_int<16>>(new_real );
                // ap_fixed<20,16> new_imag_12;
                sparse_mem_complex_t new_imag_12;
                //new_imag_12 = new_imag.range(31,16);
                new_imag_12 = new_imag.range(31,20);
                // new_imag_12.range(11,0) = new_imag.range(27,16);
                // ap_fixed<20,16> real[4] = {0};
                // ap_fixed<20,16> imag[4] = {0};
                sparse_mem_complex_t real[4] = {0};
                sparse_mem_complex_t imag[4] = {0};
                bool read_flag[4] = {0,0,0,0};
                ap_uint<9> read_coord[4] = {0};
                ap_uint<54> present_mask_ff[4] = {0};
                ap_uint<9>  min_coord_ff[4] = {0};
                switch(TX_idx){
                    case 0:
                        min_coord_ff[0] = partial_image_min_coord[1][j][k];
                        min_coord_ff[1] = partial_image_min_coord[2][j][k];
                        min_coord_ff[2] = partial_image_min_coord[3][j][k];
                        min_coord_ff[3] = partial_image_min_coord[4][j][k];
                        break;
                    case 1:
                        min_coord_ff[0] = partial_image_min_coord[0][j][k];
                        min_coord_ff[1] = partial_image_min_coord[2][j][k];
                        min_coord_ff[2] = partial_image_min_coord[3][j][k];
                        min_coord_ff[3] = partial_image_min_coord[4][j][k];
                        break;
                    case 2:
                        min_coord_ff[0] = partial_image_min_coord[0][j][k];
                        min_coord_ff[1] = partial_image_min_coord[1][j][k];
                        min_coord_ff[2] = partial_image_min_coord[3][j][k];
                        min_coord_ff[3] = partial_image_min_coord[4][j][k];
                        break;
                    case 3:
                        min_coord_ff[0] = partial_image_min_coord[0][j][k];
                        min_coord_ff[1] = partial_image_min_coord[1][j][k];
                        min_coord_ff[2] = partial_image_min_coord[2][j][k];
                        min_coord_ff[3] = partial_image_min_coord[4][j][k];
                        break;
                    case 4:
                        min_coord_ff[0] = partial_image_min_coord[0][j][k];
                        min_coord_ff[1] = partial_image_min_coord[1][j][k];
                        min_coord_ff[2] = partial_image_min_coord[2][j][k];
                        min_coord_ff[3] = partial_image_min_coord[3][j][k];
                        break;
                    default:
                        break;
                }
                // topk_update(top_real, top_imag, top_coord, new_real_12, new_imag_12, i);
                //topk_update_streaming_coord_inc(top_real, top_imag, top_coord, present_mask, valid_cnt, new_real_12, new_imag_12, i);
                topk_update_streaming_inc(top_real, top_imag, top_coord, present_mask, min_coord,valid_cnt, new_real_12, new_imag_12, i);
                // extract the 5 TX partial sums except the current TX
                // mem_read_loop1:for(int idx = 0; idx < 4; idx++){
                //     #pragma HLS UNROLL
                //     present_mask_ff[idx] = partial_image_valid_bit[other_TX_idx[idx]][j][k];
                //     read_flag[idx] = present_mask_ff[idx][i];
                //     if(read_flag[idx] == 1){
                //         //read_coord[idx] = rank_from_mask(present_mask_ff[idx], i);
                //         read_coord[idx] = rptr[idx];
                //         rptr[idx]++;
                //     }
                // }
                switch(TX_idx){
                    case 0:
                        present_mask_ff[0] = partial_image_valid_bit[1][j][k];
                        present_mask_ff[1] = partial_image_valid_bit[2][j][k];
                        present_mask_ff[2] = partial_image_valid_bit[3][j][k];
                        present_mask_ff[3] = partial_image_valid_bit[4][j][k];
                        break;
                    case 1:
                        present_mask_ff[0] = partial_image_valid_bit[0][j][k];
                        present_mask_ff[1] = partial_image_valid_bit[2][j][k];
                        present_mask_ff[2] = partial_image_valid_bit[3][j][k];
                        present_mask_ff[3] = partial_image_valid_bit[4][j][k];
                        break;
                    case 2:
                        present_mask_ff[0] = partial_image_valid_bit[0][j][k];
                        present_mask_ff[1] = partial_image_valid_bit[1][j][k];
                        present_mask_ff[2] = partial_image_valid_bit[3][j][k];
                        present_mask_ff[3] = partial_image_valid_bit[4][j][k];
                        break;
                    case 3:
                        present_mask_ff[0] = partial_image_valid_bit[0][j][k];
                        present_mask_ff[1] = partial_image_valid_bit[1][j][k];
                        present_mask_ff[2] = partial_image_valid_bit[2][j][k];   
                        present_mask_ff[3] = partial_image_valid_bit[4][j][k];
                        break;
                    case 4:
                        present_mask_ff[0] = partial_image_valid_bit[0][j][k];
                        present_mask_ff[1] = partial_image_valid_bit[1][j][k];
                        present_mask_ff[2] = partial_image_valid_bit[2][j][k];
                        present_mask_ff[3] = partial_image_valid_bit[3][j][k];
                        break;
                    default:
                        break;
                }
                for(int idx = 0; idx < 4; idx++){
                    #pragma HLS UNROLL
                    if(i >= min_coord_ff[idx] && i <= (min_coord_ff[idx] + 53)){
                        ap_uint<9> relative_k = i - min_coord_ff[idx];
                        read_flag[idx] = present_mask_ff[idx][relative_k];   
                    }
                    else{
                        read_flag[idx] = 0;
                    }
                    
                }
                for(int idx = 0; idx < 4; idx++){
                    #pragma HLS UNROLL
                    if(read_flag[idx] == 1){
                        read_coord[idx] = rptr[idx];
                        rptr[idx]++;
                    }
                }
                // for(int idx = 0; idx < 4; idx++){
                //     #pragma HLS UNROLL
                //     if(read_flag[idx] == 1){
                //         ap_int<32> temp_val = partial_image[other_TX_idx[idx]][j][k][read_coord[idx]];
                //         ap_int<16> temp_real, temp_imag;
                //         temp_real.range(15,0) = temp_val.range(31,16);
                //         temp_imag.range(15,0) = temp_val.range(15,0);
                //         real[idx].range(15,0) = temp_real.range(15, 0);
                //         imag[idx].range(15,0) = temp_imag.range(15, 0);
                //     }
                // }
                switch(TX_idx){
                    case 0:
                        real[0].range(N_bit_mem - 1,0) = partial_image[1][j][k][read_coord[0]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[0].range(N_bit_mem - 1,0) = partial_image[1][j][k][read_coord[0]].range(N_bit_mem - 1,0);
                        real[1].range(N_bit_mem - 1,0) = partial_image[2][j][k][read_coord[1]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[1].range(N_bit_mem - 1,0) = partial_image[2][j][k][read_coord[1]].range(N_bit_mem - 1,0);
                        real[2].range(N_bit_mem - 1,0) = partial_image[3][j][k][read_coord[2]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[2].range(N_bit_mem - 1,0) = partial_image[3][j][k][read_coord[2]].range(N_bit_mem - 1,0);
                        real[3].range(N_bit_mem - 1,0) = partial_image[4][j][k][read_coord[3]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[3].range(N_bit_mem - 1,0) = partial_image[4][j][k][read_coord[3]].range(N_bit_mem - 1,0);
                        break;
                    case 1:
                        real[0].range(N_bit_mem - 1,0) = partial_image[0][j][k][read_coord[0]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[0].range(N_bit_mem - 1,0) = partial_image[0][j][k][read_coord[0]].range(N_bit_mem - 1,0);
                        real[1].range(N_bit_mem - 1,0) = partial_image[2][j][k][read_coord[1]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[1].range(N_bit_mem - 1,0) = partial_image[2][j][k][read_coord[1]].range(N_bit_mem - 1,0);
                        real[2].range(N_bit_mem - 1,0) = partial_image[3][j][k][read_coord[2]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[2].range(N_bit_mem - 1,0) = partial_image[3][j][k][read_coord[2]].range(N_bit_mem - 1,0);
                        real[3].range(N_bit_mem - 1,0) = partial_image[4][j][k][read_coord[3]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[3].range(N_bit_mem - 1,0) = partial_image[4][j][k][read_coord[3]].range(N_bit_mem - 1,0);
                        break;
                    case 2:
                        real[0].range(N_bit_mem - 1,0) = partial_image[0][j][k][read_coord[0]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[0].range(N_bit_mem - 1,0) = partial_image[0][j][k][read_coord[0]].range(N_bit_mem - 1,0);
                        real[1].range(N_bit_mem - 1,0) = partial_image[1][j][k][read_coord[1]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[1].range(N_bit_mem - 1,0) = partial_image[1][j][k][read_coord[1]].range(N_bit_mem - 1,0);
                        real[2].range(N_bit_mem - 1,0) = partial_image[3][j][k][read_coord[2]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[2].range(N_bit_mem - 1,0) = partial_image[3][j][k][read_coord[2]].range(N_bit_mem - 1,0);
                        real[3].range(N_bit_mem - 1,0) = partial_image[4][j][k][read_coord[3]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[3].range(N_bit_mem - 1,0) = partial_image[4][j][k][read_coord[3]].range(N_bit_mem - 1,0);
                        break;
                    case 3:
                        real[0].range(N_bit_mem - 1,0) = partial_image[0][j][k][read_coord[0]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[0].range(N_bit_mem - 1,0) = partial_image[0][j][k][read_coord[0]].range(N_bit_mem - 1,0);
                        real[1].range(N_bit_mem - 1,0) = partial_image[1][j][k][read_coord[1]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[1].range(N_bit_mem - 1,0) = partial_image[1][j][k][read_coord[1]].range(N_bit_mem - 1,0);
                        real[2].range(N_bit_mem - 1,0) = partial_image[2][j][k][read_coord[2]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[2].range(N_bit_mem - 1,0) = partial_image[2][j][k][read_coord[2]].range(N_bit_mem - 1,0);
                        real[3].range(N_bit_mem - 1,0) = partial_image[4][j][k][read_coord[3]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[3].range(N_bit_mem - 1,0) = partial_image[4][j][k][read_coord[3]].range(N_bit_mem - 1,0);
                        break;
                    case 4:
                        real[0].range(N_bit_mem - 1,0) = partial_image[0][j][k][read_coord[0]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[0].range(N_bit_mem - 1,0) = partial_image[0][j][k][read_coord[0]].range(N_bit_mem - 1,0);
                        real[1].range(N_bit_mem - 1,0) = partial_image[1][j][k][read_coord[1]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[1].range(N_bit_mem - 1,0) = partial_image[1][j][k][read_coord[1]].range(N_bit_mem - 1,0);
                        real[2].range(N_bit_mem - 1,0) = partial_image[2][j][k][read_coord[2]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[2].range(N_bit_mem - 1,0) = partial_image[2][j][k][read_coord[2]].range(N_bit_mem - 1,0);
                        real[3].range(N_bit_mem - 1,0) = partial_image[3][j][k][read_coord[3]].range(2 * N_bit_mem - 1, N_bit_mem);
                        imag[3].range(N_bit_mem - 1,0) = partial_image[3][j][k][read_coord[3]].range(N_bit_mem - 1,0);
                        break;
                    default:
                        break;
                }
                if(read_flag[0] == 0){
                    real[0] = 0;
                    imag[0] = 0;
                }
                if(read_flag[1] == 0){
                    real[1] = 0;
                    imag[1] = 0;
                }
                if(read_flag[2] == 0){
                    real[2] = 0;
                    imag[2] = 0;
                }  
                if(read_flag[3] == 0){
                    real[3] = 0;
                    imag[3] = 0;
                }
                
                // add the new data with the other 4 TX partial sums
                ap_fixed<32, 16> final_real = new_real_12 + real[0] + real[1] + real[2] + real[3];
                ap_fixed<32, 16> final_imag = new_imag_12 + imag[0] + imag[1] + imag[2] + imag[3];
                #ifndef __SYNTHESIS__
                    if(TX_idx == 4){
                        if(final_real != 0 || final_imag != 0){
                            debug_vector_adder << "At line " << j << ", phi " << k << ", point " << i << ": New real = " << new_real_12 << ", New imag = " << new_imag_12
                                        << " | Other TX reals = [" << real[0] << ", " << real[1] << ", " << real[2] << ", " << real[3] << "]"
                                        << " | Other TX imags = [" << imag[0] << ", " << imag[1] << ", " << imag[2] << ", " << imag[3] << "]"
                                        << " | Final real = " << final_real << ", Final imag = " << final_imag << std::endl;
                
                        }
                    }
                #endif
                bf_data_t_64 output;
                output = final_real * final_real + final_imag * final_imag;
                out_data.write(output);
                SA_partial_image.write(new_val); // pass through the SA_partial_image stream
                // debug_vector_adder << "Final output at line " << j << ", phi " << k << ", point " << i << ": " << output << std::endl;
            }
             // update the partial_image with the new top K values
            partial_image_valid_bit[TX_idx][j][k] = present_mask;
            partial_image_min_coord[TX_idx][j][k] = min_coord;
             for(int p = 0; p < K; p++){
                // partial_image[TX_idx][j][k][p].range(31,16) = top_real[p].range(15, 0);
                // partial_image[TX_idx][j][k][p].range(15,0)  = top_imag[p].range(15, 0);
                partial_image[TX_idx][j][k][p].range(23,12) = top_real[p].range(11, 0);
                partial_image[TX_idx][j][k][p].range(11,0)  = top_imag[p].range(11, 0);
                // partial_image_coord[TX_idx][j][k][p] = top_coord[p];
                // output the top K values to out_data stream
                #ifndef __SYNTHESIS__
                    debug_topk << "TX " << (int)(TX_idx) << ", Line " << j << ", Phi " << k << ", Top " << p << ": Real = " << top_real[p] << ", Imag = " << top_imag[p] << ", Coord = " << top_coord[p] << " \n ";
                    debug_topk << "min coord: " << min_coord << ", present_mask: " << present_mask.to_string(2) << "\n";
                #endif
            }
        }
    }
    //debug_vector_adder.close();
     #ifndef __SYNTHESIS__
        std::cout << "vector_adder processing completed." << std::endl;
        debug_vector_adder.close();
        debug_topk.close();
    #endif
}
//////////////////////////////////////
//   sparce memory access pattern   //
//////////////////////////////////////
template<int num_bit>
bool topk_update(ap_int<num_bit> top_real[K], ap_int<num_bit> top_imag[K], ap_uint<9> top_coord[K], ap_int<num_bit> x_real, ap_int<num_bit> x_imag, ap_uint<9> x_coord){
	#pragma HLS INLINE
    
    ap_int<2*num_bit + 1> x_mag_sq = x_real * x_real + x_imag * x_imag;
    ap_int<2*num_bit + 1> top_mag_least = top_real[K-1] * top_real[K-1] + top_imag[K-1] * top_imag[K-1];
    // 如果新值不大於目前第 K 大，直接返回
    if (x_mag_sq <= top_mag_least) return false;

    // 插入排序式更新：由尾往前推
    int pos = K - 1;
    for (int i = K - 1; i > 0; --i) {
        #pragma HLS UNROLL
        if (x_mag_sq > (top_real[i-1] * top_real[i-1] + top_imag[i-1] * top_imag[i-1])) {
            top_real[i] = top_real[i-1];
            top_imag[i] = top_imag[i-1];
            top_coord[i] = top_coord[i-1];
            pos = i - 1;
        } else {
            break;
        }
    }
    top_real[pos] = x_real;
    top_imag[pos] = x_imag;
    top_coord[pos] = x_coord;

    return true;
}
// topK update no control
void topk_update_streaming_inc(
    sparse_mem_complex_t   top_real[K],
    sparse_mem_complex_t   top_imag[K],
    ap_uint<9>   top_coord[K],      // 有效資料永遠 compact 在前面，且 coord 遞增
    ap_uint<54>   &present_mask,
    ap_uint<9>   &min_coord,
    ap_uint<6>   &valid_cnt,        // 0..K
    sparse_mem_complex_t   x_real,
    sparse_mem_complex_t   x_imag,
    ap_uint<9>   x_coord            // 0..283，且每次呼叫都遞增
){
    ap_int<33> x_mag = x_real*x_real + x_imag*x_imag;
    ap_int<33> x_mag_init = top_real[0]*top_real[0] + top_imag[0]*top_imag[0];
    ap_uint<9> min_next_coord = min_coord;
    for(int i = 1; i < 53; i++){ // find the next min coord in case we need to evict the current min coord, the max distance between current min coord and next min coord is 53 because of the topk update rule
        #pragma HLS UNROLL
        if(present_mask[i] == 1){
            min_next_coord = min_coord + i;
            break;
        }
    }
    // if(valid_cnt < K && x_mag != 0){
    //     top_real[valid_cnt] = x_real;
    //     top_imag[valid_cnt] = x_imag;
    //     top_coord[valid_cnt] = x_coord;
    //     present_mask[x_coord] = 1;
    //     valid_cnt++;
    //     return;
    // }
    if(valid_cnt == 0 && x_mag != 0){
        top_real[0] = x_real;
        top_imag[0] = x_imag;
        top_coord[0] = x_coord;
        present_mask[0] = 1;
        min_coord = x_coord;
        valid_cnt++;
        return;
    }
    else if(x_coord <= (min_coord + 53) && valid_cnt < K && x_mag != 0){
        top_real[valid_cnt] = x_real;
        top_imag[valid_cnt] = x_imag;
        top_coord[valid_cnt] = x_coord;
        present_mask[x_coord - min_coord] = 1;
        valid_cnt++;
        return;
    }
    else if(x_coord <= (min_coord + 53) && valid_cnt == (K - 1) && x_mag > x_mag_init){
        for(int i = 0; i < K - 1; i++){
            #pragma HLS UNROLL
            top_real[i] = top_real[i+1];
            top_imag[i] = top_imag[i+1];
            top_coord[i] = top_coord[i+1];
        }
        top_real[K - 1] = x_real;
        top_imag[K - 1] = x_imag;
        top_coord[K - 1] = x_coord;
        present_mask[x_coord - min_coord] = 1;
        min_coord = min_next_coord;
        return;
    }
}
// topk update but use coord to sort and keep gate with mask
void topk_update_streaming_coord_inc(
    ap_int<16>   top_real[K],
    ap_int<16>   top_imag[K],
    ap_uint<9>   top_coord[K],      // 有效資料永遠 compact 在前面，且 coord 遞增
    ap_uint<284> &present_mask,
    ap_uint<6>   &valid_cnt,        // 0..K
    ap_int<16>   x_real,
    ap_int<16>   x_imag,
    ap_uint<9>   x_coord            // 0..283，且每次呼叫都遞增
){
#pragma HLS INLINE off

    // 不重複（理論上不會發生）
    

    ap_int<33> x_mag = x_real*x_real + x_imag*x_imag;
    if (present_mask[x_coord] || x_mag == 0) return;
    // 1) 還沒滿：直接 append 到尾巴（因為 coord 遞增）
    if (valid_cnt < K) {
        top_real[valid_cnt]  = x_real;
        top_imag[valid_cnt]  = x_imag;
        top_coord[valid_cnt] = x_coord;
        present_mask[x_coord] = 1;
        valid_cnt++;
        #ifndef __SYNTHESIS__
            std::cout << "Inserted point at coord " << x_coord << " into position " << valid_cnt - 1 << " (not full yet)." << std::endl;
        #endif
        return;
    }

    // 2) 已滿：找最小 magnitude 那筆
    int min_pos = 0;
    ap_int<33> min_mag = top_real[0]*top_real[0] + top_imag[0]*top_imag[0];
    for (int i = 1; i < K; i++) {
        #pragma HLS UNROLL
        ap_int<33> m = top_real[i]*top_real[i] + top_imag[i]*top_imag[i];
        if (m < min_mag) { min_mag = m; min_pos = i; }
    }

    // Gate：比最小還小就丟掉
    if (x_mag <= min_mag ) return;

    // 3) 踢掉 min_pos：mask 清掉 + 左移補洞（保持 compact & coord 遞增）
    ap_uint<9> evict_coord = top_coord[min_pos];
    present_mask[evict_coord] = 0;

    for (int i = 0; i < K-1; i++) {
        #pragma HLS UNROLL
        if (i >= min_pos) {
            top_real[i]  = top_real[i+1];
            top_imag[i]  = top_imag[i+1];
            top_coord[i] = top_coord[i+1];
        }
    }

    // 4) 新點放最後（coord 最大）
    top_real[K-1]  = x_real;
    top_imag[K-1]  = x_imag;
    top_coord[K-1] = x_coord;
    present_mask[x_coord] = 1;
}

ap_uint<9> rank_from_mask(const ap_uint<284> mask, ap_uint<9> coord) {
    #pragma HLS inline

    ap_uint<9> count = 0;
    for (int i = 0; i < 248; i++){
        #pragma HLS UNROLL
        if (i < coord){
            count += mask[i];
        }
    }
    return count;
}