// Eoin O'Connell
// Feb 11 2025
// eoconnell@hmc.edu

// NEED TO FIX:

// need to fix sign for overflow (is the infinity part of the add or mul?)
// underflow for subnormal
// number of leading zeros to determine subnormal


module fma16 (
    input  logic [15:0] x, y, z,
    input  logic        mul, add, negp, negz,
    input  logic [1:0]  roundmode,
    output logic [15:0] result,
    output logic [3:0]  flags
);
    
    // Handle input configurations (mul, add, negp, negz)
    logic [15:0] y_input, z_input;
    assign y_input = mul ? {negp ^ y[15], y[14:0]} : 16'h3c00;
    assign z_input = add ? {negz ^ z[15], z[14:0]} : 16'b0;
    


    // Initiallizing internal logic signals

    // for fma algorithm
    logic [4:0] X_e, Y_e, Z_e;
    logic [5:0] P_e, M_e, shift_amount;
    logic X_nonzero, Y_nonzero, Z_nonzero, X_subnorm, Y_subnorm, Z_subnorm, small_res, subtract, sign, P_sign;
    logic [10:0] X_m, Y_m, Z_m;
    logic [9:0] M_frac;
    logic [21:0] P_m, P_intermediate, A_m, Z_m_padded, P_m_shifted;
    logic [22:0] S_m, intermediate_m;
    logic [4:0] P_zeros, P_shift_amount;
    logic [4:0] leading_zeros;
    logic [4:0] P_shift;
    logic [40:0] S_shift, res_shift;
    logic [15:0] fma_result;  
    logic res_subnorm, smaller, negligable;
    logic small_product;

    // for flags and special cases
    logic invalid, overflow, underflow, inexact, overflow_maxnum;
    logic sNaN, qNaN, zero_times_infinity, infinity_minus_infinity;
    logic X_max_exp, Y_max_exp, Z_max_exp;
    logic X_inf, Y_inf, Z_inf;
    logic valid_inf, inf_sign;
    logic exact_zero, product_nonzero;

    // for rounding
    logic [15:0] rounded;
    logic round_overflow;

        
    // SIGN LOGIC
    assign P_sign = x[15] ^ y_input[15];  // sign is 1 when product is negative
    assign subtract = P_sign ^ z_input[15];     // subtract is 0 when product and z have same sign


    // Process Inputs and detect special cases
    assign X_e = x[14:10];
    assign Y_e = y_input[14:10];
    assign Z_e = z_input[14:10];

    assign X_nonzero = |x[14:0];
    assign Y_nonzero = |y_input[14:0];
    assign Z_nonzero = |z_input[14:0];

    assign X_inf = (X_e == 31) && (x[9:0] == 10'b0);
    assign Y_inf = (Y_e == 31) && (y_input[9:0] == 10'b0);
    assign Z_inf = (Z_e == 31) && (z_input[9:0] == 10'b0);

    assign X_subnorm = (X_e == 0) & (X_nonzero);
    assign Y_subnorm = (Y_e == 0) & (Y_nonzero);
    assign Z_subnorm = (Z_e == 0) & (Z_nonzero);

    assign X_m = {X_nonzero & (~X_subnorm), x[9:0]};
    assign Y_m = {Y_nonzero & (~Y_subnorm), y_input[9:0]};
    assign Z_m = {Z_nonzero & (~Z_subnorm), z_input[9:0]};

    assign X_max_exp = X_e == 31;
    assign Y_max_exp = Y_e == 31;
    assign Z_max_exp = Z_e == 31;

    assign sNaN = (X_max_exp && !x[9] && |x[8:0]) || 
        (Y_max_exp && !y_input[9] && |y_input[8:0]) || 
        (Z_max_exp && !z_input[9] && |z_input[8:0]);
    assign qNaN = (X_max_exp && x[9] && |x[9:0]) || 
        (Y_max_exp && y_input[9] && |y_input[9:0]) || 
        (Z_max_exp && z_input[9] && |z_input[9:0]);        
    assign zero_times_infinity = (X_inf & !Y_nonzero) | (!X_nonzero & Y_inf);
    assign infinity_minus_infinity = (P_sign ^ z[15]) & (X_inf | Y_inf) & Z_inf;
    assign invalid = sNaN || zero_times_infinity || infinity_minus_infinity;

    assign valid_inf = (X_inf || Y_inf || Z_inf) && !(zero_times_infinity || infinity_minus_infinity);
    assign inf_sign = ((X_inf && P_sign) || (Y_inf && P_sign) || (Z_inf && z[15]));

    assign product_nonzero = (X_nonzero && Y_nonzero);

    //$display("X: %b, Y: %b, X_m: %b, Y_m: %b", x, y, X_m, Y_m);




    

    // Starting FMA algorithm
    assign P_m = X_m * Y_m;
    // FIX SMALL PRODUCT
    assign small_product = (X_e + Y_e + {5'b0, P_m[21]} + {5'b0, X_subnorm} + {5'b0, Y_subnorm}) < 15;
    assign P_e = (product_nonzero ? (X_e + Y_e - 15 + {5'b0, X_subnorm} + {5'b0, Y_subnorm} + {5'b0, P_m[21]}) : 0);
    assign P_m_shifted = P_m[21] ? P_m : P_m << 1;    

    assign Z_m_padded = {Z_m, 11'b0};

    always_comb begin
        // smaller = 1 when product greater than addend
        if (P_e == {1'b0, Z_e}) begin
            smaller = !small_product && (P_m_shifted > Z_m_padded);
        end else begin
            smaller = !small_product && (P_e > {1'b0, Z_e});
        end

        // determine shift amount for smaller component
        if (smaller) begin
            shift_amount = (P_e - (Z_e + {5'b0, Z_subnorm})); 
        end else begin
            shift_amount = (Z_e + {5'b0, Z_subnorm}) - P_e; 
        end
    end
        
    // NEED TO FIX P_shift AND USE
    //assign P_shift = small_product ? (Z_e + {4'b0, Z_subnorm}) + (X_e + Y_e + {4'b0, X_subnorm} + {4'b0, Y_subnorm}) : (Z_e + {4'b0, Z_subnorm}) - P_e;
    assign A_m = smaller ? Z_m_padded >> shift_amount : P_m_shifted >> shift_amount; 

    logic [4:0] count_one, count_two;
    logic [21:0] not_shifted, shifted;

    always_comb begin
        if (smaller) begin
            not_shifted = P_m_shifted;
            shifted = Z_m_padded;
        end else begin
            not_shifted = Z_m_padded;
            shifted = P_m_shifted;
        end
    end

    count_ones c1(A_m, count_one);
    count_ones c2(shifted, count_two);
    assign negligable = ~(count_one == count_two) && !A_m[0];
    //assign negligable = smaller ? (Z_m_padded != 0) && (A_m == 0) : (P_m_shifted != 0) && (A_m == 0);

    always_comb begin
        if (subtract) begin
            S_m = smaller ? not_shifted - A_m - {22'b0, negligable} : not_shifted - A_m - {22'b0, negligable}; 
        end else begin
            S_m = A_m + not_shifted + {22'b0, negligable};
        end
    end

    lzd23 l(S_m, leading_zeros);
    

        

    // Check for exact zero results (not rounded to zero)
    assign exact_zero = !(|S_m);




    // NEED LZD DETECTOR:
    // need to check for subnormal numbers: somthing like (leading_zeros > 31) res_e = small_res ? 5'd0 : 5'(30 - leading_zeros);

    // need to add two because I added extra bit for multiplication, and then another extra bit for addition 
    
    assign res_subnorm = leading_zeros > 13;  // FIX THIS NUMBER, RANDOM GUESS, NEED TO INCORPORATE EXPONENT AS WELL
    
    assign M_e = res_subnorm ? 0 : (smaller ? P_e + 1 - leading_zeros: Z_e + 1 - leading_zeros);
    assign overflow = M_e > 30;
    assign intermediate_m = S_m << (res_subnorm ? 13 : leading_zeros);  // FIX THIS NUMBER (13), RANDOM GUESS
    assign M_frac = intermediate_m[21:12];

    assign sign = ((smaller) ? P_sign : z_input[15]) && !exact_zero || (exact_zero && ((roundmode == 2'b10 && subtract) || (!subtract && z[15]))); // roundmode == 2'b10

    assign fma_result = {sign, M_e[4:0], M_frac};  // deleted   (sign & (~small_res))
    



    assign underflow = (M_e == 0 && |intermediate_m[11:0]); 



    // initialize rounding module
    rounding r(sign, M_e[4:0], intermediate_m, roundmode, rounded, round_overflow);



    
    
        
    // Flags
    assign inexact = overflow || |intermediate_m[11:0] || (negligable && subtract); // find the range of bits that is after mantissa, (or of all bits, if any is 1, then inexact)

    // flags and special cases
    assign flags = {invalid, (overflow || round_overflow) && !invalid && !qNaN && !valid_inf, underflow && !invalid && !qNaN && !valid_inf, (inexact || negligable) && !invalid && !qNaN && !valid_inf};

    // check if the rounding mode + sign combination makes overflow be maxnum instead of infinity
    assign overflow_maxnum = sign ? roundmode == 2'b00 || roundmode == 2'b11 : roundmode == 2'b00 || roundmode == 2'b10;
    
    always_comb begin

        // set final result
        casez(flags)
            4'b1???: result = 16'b0_11111_1000000000;
            4'b010?: result = {sign, overflow_maxnum ? 15'b11110_1111111111 : 15'b11111_0000000000};  // depends on rounding mode
            default: result = qNaN ? 16'b0_11111_1000000000 : (valid_inf ? {inf_sign, 15'b11111_0000000000} : rounded);
        endcase
        
        if (x == 16'h0401 && y == 16'h7801 && z == 16'h7bff) begin
            //$display("intermediate %b", intermediate_m);
            //$display("x: %b, y: %b, z: %b", x, y_input, z_input);
            //$display("P_m: %b, P_e: %b, A_m: %b, Z_m_padded: %b, S_m: %b, intermediate: %b, smaller: %b, negligable: %b", P_m, P_e, A_m, Z_m_padded, S_m, intermediate_m, smaller, negligable);
            //$display("%b", {1'b0, Z_m_padded});
            //$display("%b", {1'b0, A_m});
            //$display("%b", S_m);
            //$display("overflow reg: %b, overflow round: %b", overflow, round_overflow);
            //$display("count one: %b, count two: %b, not shifted: %b, overflow_maxnum: %b", count_one, count_two, not_shifted, overflow_maxnum);
            //$display("P_e: %b, Z_e: %b, res: %b, shift amount: %b", P_e, Z_e, (P_e - (Z_e + {5'b0, Z_subnorm})), shift_amount);
        end
        
    end

endmodule





module rounding(
    input logic sign, 
    input logic [4:0] exponent, 
    input logic [22:0] intermediate_m, 
    input logic [1:0] roundmode,
    output logic [15:0] result,
    output logic overflow);
    
    logic [9:0] mantissa, final_mantissa, truncated; // do I want leading 1?
    logic [9:0] rounded_up;
    logic increase_exp;
    logic [9:0] rounded_down;
    logic [4:0] round_up_exp, final_exponent;
    logic [11:0] extra_bits;
    logic odd, rne_up;
    logic inexact;
    logic rounded_sign;
    logic round_up;

    /*
    Notes:
    00: round to zero    (truncate)
    01: round to even       
    10: round down (toward negative infinity)
    11: round up (toward positive infinity)

    Overflow:
    */

    assign mantissa = intermediate_m[21:12];

    assign increase_exp = mantissa == 10'b1111111111;
    assign rounded_up = increase_exp ? 10'b0000000000 : mantissa + {10'b0000000001};
    assign round_up_exp = exponent + {4'b0000, increase_exp};
    assign truncated = mantissa;
    assign extra_bits = intermediate_m[11:0];

    assign inexact = |extra_bits[11:0];

    // is this the correct check for even/odd?
    assign odd = mantissa[0];
    assign rne_up = extra_bits[11] && (|extra_bits[10:0] || odd);  

    always_comb begin
        /*
        case (roundmode)
            2'b00:  begin 
                        final_mantissa = truncated;
                        final_exponent = exponent;
                        end
            2'b01:  begin 
                        final_mantissa = rne_up ? rounded_up : truncated;
                        final_exponent = rne_up ? round_up_exp : exponent;
                        end
            2'b10:  begin
                        final_mantissa = sign && inexact ? rounded_up : truncated;
                        final_exponent = sign && inexact ? round_up_exp : exponent;
                        end
            default: begin
                        final_mantissa = !sign && inexact ? rounded_up : truncated;
                        final_exponent = !sign && inexact ? round_up_exp : exponent;
                        end
        endcase
        */
        case (roundmode)
            2'b00:  begin 
                        round_up = 0;
                        end
            2'b01:  begin 
                        round_up = rne_up;
                        end
            2'b10:  begin
                        round_up = sign && inexact;
                        end
            default: begin
                        round_up = !sign && inexact;
                        end
        endcase
        final_mantissa = round_up ? rounded_up : truncated;
        final_exponent = round_up ? round_up_exp : exponent;
        if (result == 16'h7801) begin
            //$display("intermediate_m %b, truncated %b, roundedup %b, result %b", intermediate_m, truncated, rounded_up, final_mantissa);
        end
    end
        
    assign overflow = ((exponent == 5'b11110) && increase_exp && round_up) || (exponent == 5'b11111);

    assign result = {sign, final_exponent, final_mantissa};    

endmodule





module lzd23(input logic [22:0] a,
            output logic [4:0] leading_zeros);

    always_comb begin
        casez (a)
            23'b00000000000000000000000: leading_zeros = 5'd23;
            23'b00000000000000000000001: leading_zeros = 5'd22;
            23'b0000000000000000000001?: leading_zeros = 5'd21;
            23'b000000000000000000001??: leading_zeros = 5'd20;
            23'b00000000000000000001???: leading_zeros = 5'd19;
            23'b0000000000000000001????: leading_zeros = 5'd18;
            23'b000000000000000001?????: leading_zeros = 5'd17;
            23'b00000000000000001??????: leading_zeros = 5'd16;
            23'b0000000000000001???????: leading_zeros = 5'd15;
            23'b000000000000001????????: leading_zeros = 5'd14;
            23'b00000000000001?????????: leading_zeros = 5'd13;
            23'b0000000000001??????????: leading_zeros = 5'd12;
            23'b000000000001???????????: leading_zeros = 5'd11;
            23'b00000000001????????????: leading_zeros = 5'd10;
            23'b0000000001?????????????: leading_zeros = 5'd9;
            23'b000000001??????????????: leading_zeros = 5'd8;
            23'b00000001???????????????: leading_zeros = 5'd7;
            23'b0000001????????????????: leading_zeros = 5'd6;
            23'b000001?????????????????: leading_zeros = 5'd5;
            23'b00001??????????????????: leading_zeros = 5'd4;
            23'b0001???????????????????: leading_zeros = 5'd3;
            23'b001????????????????????: leading_zeros = 5'd2;
            23'b01?????????????????????: leading_zeros = 5'd1;
            default:                    leading_zeros = 5'd0;
        endcase

    end
endmodule

module count_ones (
    input  logic [21:0] in,
    output logic [4:0]  count
);
    assign count =
        in[0]  + in[1]  + in[2]  + in[3]  + in[4]  + in[5]  + in[6]  + in[7]  +
        in[8]  + in[9]  + in[10] + in[11] + in[12] + in[13] + in[14] + in[15] +
        in[16] + in[17] + in[18] + in[19] + in[20] + in[21];
endmodule


