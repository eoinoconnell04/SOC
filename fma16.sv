// Eoin O'Connell
// Feb 11 2025
// eoconnell@hmc.edu

// NEED TO FIX:

// need to fix sign for overflow (is the infinity part of the add or mul?)
// underflow for subnormal
// number of leading zeros to determine subnormal


// Timing: 7.727844
// Area: 29180.434504


module fma16 (
    input  logic [15:0] x, y, z,
    input  logic        mul, add, negp, negz,
    input  logic [1:0]  roundmode,
    output logic [15:0] result,
    output logic [3:0]  flags
);
    
    // Handle input configurations (mul, add, negp, negz)
    logic [15:0] Y_input, Z_input;
    operation op(mul, add, negp, negz, y, z, Y_input, Z_input);
    


    // Initiallizing internal logic signals

    // for fma algorithm
    logic [4:0] X_e, Y_e, Z_e;
    logic [5:0] P_e, M_e, shift_amount;
    logic X_nonzero, Y_nonzero, Z_nonzero, X_subnorm, Y_subnorm, Z_subnorm, subtract, sign, P_sign;
    logic [10:0] X_m, Y_m, Z_m;
    logic [9:0] M_frac;
    logic [21:0] P_m, P_intermediate, A_m, Z_m_padded, P_m_shifted;
    logic [22:0] S_m, intermediate_m;
    logic [4:0] P_zeros, P_shift_amount;
    logic [4:0] leading_zeros;
    logic res_subnorm, smaller, negligable;
    logic small_product;

    // for flags and special cases
    logic invalid, overflow, underflow, inexact, overflow_maxnum, qNaN;
    logic valid_inf, inf_sign;
    logic exact_zero, product_nonzero;

    // for rounding
    logic [15:0] rounded;
    logic round_overflow;

        
    // SIGN LOGIC
    assign P_sign = x[15] ^ Y_input[15];  // sign is 1 when product is negative
    assign subtract = P_sign ^ Z_input[15];     // subtract is 0 when product and z have same sign


    // Process Inputs to be used for special cases and FMA algorithm
    inputs i(x, Y_input, Z_input, X_e, Y_e, Z_e, X_nonzero, Y_nonzero, X_subnorm, Y_subnorm, Z_subnorm, X_m, Y_m, Z_m); 
    /*
    assign X_e = x[14:10];
    assign Y_e = Y_input[14:10];
    assign Z_e = Z_input[14:10];

    assign X_nonzero = |x[14:0];
    assign Y_nonzero = |Y_input[14:0];
    assign Z_nonzero = |Z_input[14:0];

    assign X_subnorm = (X_e == 0) & (X_nonzero);
    assign Y_subnorm = (Y_e == 0) & (Y_nonzero);
    assign Z_subnorm = (Z_e == 0) & (Z_nonzero);

    assign X_m = {X_nonzero & (~X_subnorm), x[9:0]};
    assign Y_m = {Y_nonzero & (~Y_subnorm), Y_input[9:0]};
    assign Z_m = {Z_nonzero & (~Z_subnorm), Z_input[9:0]};
    */

    special_cases sp(x, Y_input, Z_input, X_e, Y_e, Z_e, X_nonzero, Y_nonzero, P_sign, valid_inf, inf_sign, invalid, qNaN);

    assign product_nonzero = (X_nonzero && Y_nonzero);


    // Starting FMA algorithm

    // Find Product P_m
    product p(X_m, Y_m, X_e, Y_e, X_subnorm, Y_subnorm, product_nonzero, small_product, P_e, P_m_shifted);

    /*
    assign P_m = X_m * Y_m;
    assign small_product = (X_e + Y_e + {5'b0, P_m[21]} + {5'b0, X_subnorm} + {5'b0, Y_subnorm}) < 15;
    assign P_e = (product_nonzero ? (X_e + Y_e - 15 + {5'b0, X_subnorm} + {5'b0, Y_subnorm} + {5'b0, P_m[21]}) : 0);
    assign P_m_shifted = P_m[21] ? P_m : P_m << 1;    
    */

    logic [21:0] not_shifted, shifted;
    align_shift s(P_m_shifted, P_e, Z_m, Z_e, small_product, Z_subnorm, smaller, shift_amount, not_shifted, shifted, A_m);
    //always_comb begin $display("%b, %b, %b, %b, %b, %b, %b, %b, %b, %b, %b", P_m_shifted, P_e, Z_m, Z_e, small_product, Z_subnorm, smaller, shift_amount, not_shifted, shifted, A_m); end
    /*
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
        
    assign A_m = smaller ? Z_m_padded >> shift_amount : P_m_shifted >> shift_amount; 

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
    */

    logic prec_lost;
    precision_lost pl(shifted, shift_amount, prec_lost);

    // negligable occurs when there is a loss of precision that could chagne the last bit of A_m from a 0 to a 1
    assign negligable = prec_lost /*~(count_one == count_two)*/ && !A_m[0];

    add_sub as(A_m, not_shifted, smaller, negligable, subtract, S_m);

    lzd23 l(S_m, leading_zeros);
       

    // Check for exact zero results (not rounded to zero)
    assign exact_zero = !(|S_m);

    // need to check for subnormal numbers: somthing like (leading_zeros > 31) res_e = small_res ? 5'd0 : 5'(30 - leading_zeros);    
    assign res_subnorm = leading_zeros > 13;  // FIX THIS NUMBER, RANDOM GUESS, NEED TO INCORPORATE EXPONENT AS WELL
    
    assign M_e = res_subnorm ? 0 : (smaller ? P_e + 1 - leading_zeros: Z_e + 1 - leading_zeros);
    assign overflow = M_e > 30;
    assign intermediate_m = S_m << (res_subnorm ? 13 : leading_zeros); 
    assign M_frac = intermediate_m[21:12];

    assign sign = ((smaller) ? P_sign : Z_input[15]) && !exact_zero || (exact_zero && ((roundmode == 2'b10 && subtract) || (!subtract && Z_input[15]))); // roundmode == 2'b10    



    assign underflow = (M_e == 0 && |intermediate_m[11:0]); 

    // initialize rounding module
    rounding r(sign, M_e[4:0], intermediate_m, roundmode, rounded, round_overflow);

    // Flags
    assign inexact = overflow || |intermediate_m[11:0] || (negligable && subtract); // find the range of bits that is after mantissa, (or of all bits, if any is 1, then inexact)

    // flags and special cases
    assign flags = {invalid, (overflow || round_overflow) && !invalid && !qNaN && !valid_inf, underflow && !invalid && !qNaN && !valid_inf, (inexact || negligable) && !invalid && !qNaN && !valid_inf};

    // check if the rounding mode + sign combination makes overflow be maxnum instead of infinity
    assign overflow_maxnum = sign ? roundmode == 2'b00 || roundmode == 2'b11 : roundmode == 2'b00 || roundmode == 2'b10;
    
    final_result fr(sign, overflow_maxnum, qNaN, valid_inf, inf_sign, rounded, flags, result);
    
endmodule

module final_result(
    input logic sign, overflow_maxnum, qNaN, valid_inf, inf_sign, 
    input logic [15:0] rounded,
    input logic [3:0] flags,
    output logic [15:0] result
);
    always_comb begin

        // set final result
        casez(flags)
            4'b1???: result = 16'b0_11111_1000000000;
            4'b010?: result = {sign, overflow_maxnum ? 15'b11110_1111111111 : 15'b11111_0000000000};  // depends on rounding mode
            default: result = qNaN ? 16'b0_11111_1000000000 : (valid_inf ? {inf_sign, 15'b11111_0000000000} : rounded);
        endcase
        
    end
endmodule

// performs rounding given an unrounded answer combined with a specified rounding mode
/*
    Rounding Modes:
    00: round to zero    (truncate)
    01: round to even       
    10: round down (toward negative infinity)
    11: round up (toward positive infinity)
*/
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

    end
        
    assign overflow = ((exponent == 5'b11110) && increase_exp && round_up) || (exponent == 5'b11111);

    assign result = {sign, final_exponent, final_mantissa};    

endmodule




// leading zeros detector
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

// selects and modifies inputs based on the operational inputs
module operation (
    input logic mul, add, negp, negz, 
    input logic [15:0] y, z, 
    output logic [15:0] Y_input, Z_input);
    assign Y_input = mul ? {negp ^ y[15], y[14:0]} : 16'h3c00;
    assign Z_input = add ? {negz ^ z[15], z[14:0]} : 16'b0;
endmodule

module product(
    input logic [10:0] X_m, Y_m, 
    input logic [4:0] X_e, Y_e, 
    input logic X_subnorm, Y_subnorm, product_nonzero,
    output logic small_product,
    output logic [5:0] P_e,
    output logic [21:0] P_m_shifted);
    logic [21:0] P_m;
    assign P_m = X_m * Y_m;
    assign small_product = (X_e + Y_e + {5'b0, P_m[21]} + {5'b0, X_subnorm} + {5'b0, Y_subnorm}) < 15;
    assign P_e = (product_nonzero ? (X_e + Y_e - 15 + {5'b0, X_subnorm} + {5'b0, Y_subnorm} + {5'b0, P_m[21]}) : 0);
    assign P_m_shifted = P_m[21] ? P_m : P_m << 1;  
endmodule


// indentifies special cases
module special_cases(
    input logic [15:0] x, Y_input, Z_input, 
    input logic [4:0] X_e, Y_e, Z_e, 
    input logic X_nonzero, Y_nonzero, P_sign,
    output logic valid_inf, inf_sign, invalid, qNaN);

    logic sNaN, zero_times_infinity, infinity_minus_infinity;
    logic X_or_8bits, Y_or_8bits, Z_or_8bits;
    logic X_max_exp, Y_max_exp, Z_max_exp;
    logic X_inf, Y_inf, Z_inf;

    assign X_max_exp = X_e == 31;
    assign Y_max_exp = Y_e == 31;
    assign Z_max_exp = Z_e == 31;
    assign X_or_8bits = |x[8:0];
    assign Y_or_8bits = |Y_input[8:0];
    assign Z_or_8bits = |Z_input[8:0];

    assign X_inf = (X_max_exp) && (x[9:0] == 10'b0);
    assign Y_inf = (Y_max_exp) && (Y_input[9:0] == 10'b0);
    assign Z_inf = (Z_max_exp) && (Z_input[9:0] == 10'b0);

    assign sNaN = (X_max_exp && !x[9] && X_or_8bits) || 
        (Y_max_exp && !Y_input[9] && Y_or_8bits) || 
        (Z_max_exp && !Z_input[9] && Z_or_8bits);
    assign qNaN = (X_max_exp && x[9] && (x[9] || X_or_8bits)) || 
        (Y_max_exp && Y_input[9] && (Y_input[9] || Y_or_8bits)) || 
        (Z_max_exp && Z_input[9] && (Z_input[9] || Z_or_8bits));        
    assign zero_times_infinity = (X_inf & !Y_nonzero) | (!X_nonzero & Y_inf);
    assign infinity_minus_infinity = (P_sign ^ Z_input[15]) & (X_inf | Y_inf) & Z_inf;
    assign invalid = sNaN || zero_times_infinity || infinity_minus_infinity;

    assign valid_inf = (X_inf || Y_inf || Z_inf) && !(zero_times_infinity || infinity_minus_infinity);
    assign inf_sign = ((X_inf && P_sign) || (Y_inf && P_sign) || (Z_inf && Z_input[15]));

endmodule

// performs the addition and subtraction operations
module add_sub(
    input logic [21:0] A_m, not_shifted, 
    input logic smaller, negligable, subtract, 
    output logic [22:0] S_m);
    always_comb begin
        if (subtract) begin
            S_m = smaller ? not_shifted - A_m - {22'b0, negligable} : not_shifted - A_m - {22'b0, negligable}; 
        end else begin
            S_m = A_m + not_shifted + {22'b0, negligable};
        end
    end
endmodule

// returns a 1 whenever the right shift of the input by the shift amount would shift away a 1
module precision_lost(
    input  logic [21:0] in,
    input  logic [5:0] shift_amount,
    output logic prec_lost);
    logic [21:0] shifted;
    logic far_shift;
    assign far_shift = shift_amount > 6'd22;
    assign shifted = in << (far_shift ? 0 : 6'd22-shift_amount);
    assign prec_lost = |shifted;
    //always_comb begin $display("in %b, shift_amount: %d, shifted %b, far_shift %b, prec_lost %b", in, shift_amount, shifted, far_shift, prec_lost); end
endmodule

module inputs(
    input logic [15:0] x, Y_input, Z_input, 
    output logic [4:0] X_e, Y_e, Z_e, 
    output logic X_nonzero, Y_nonzero, X_subnorm, Y_subnorm, Z_subnorm, 
    output logic [10:0] X_m, Y_m, Z_m); 
    logic Z_nonzero;
    assign X_e = x[14:10];
    assign Y_e = Y_input[14:10];
    assign Z_e = Z_input[14:10];

    assign X_nonzero = |x[14:0];
    assign Y_nonzero = |Y_input[14:0];
    assign Z_nonzero = |Z_input[14:0];

    assign X_subnorm = (X_e == 0) & (X_nonzero);
    assign Y_subnorm = (Y_e == 0) & (Y_nonzero);
    assign Z_subnorm = (Z_e == 0) & (Z_nonzero);

    assign X_m = {X_nonzero & (~X_subnorm), x[9:0]};
    assign Y_m = {Y_nonzero & (~Y_subnorm), Y_input[9:0]};
    assign Z_m = {Z_nonzero & (~Z_subnorm), Z_input[9:0]};
endmodule

module align_shift(
    input logic [21:0] P_m_shifted,
    input logic [5:0] P_e,
    input logic [10:0] Z_m,
    input logic [4:0] Z_e,
    input logic small_product, Z_subnorm,
    output logic smaller,
    output logic [5:0] shift_amount,
    output logic [21:0] not_shifted, shifted, A_m
    );
    logic [21:0] Z_m_padded;
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
        
    assign A_m = smaller ? Z_m_padded >> shift_amount : P_m_shifted >> shift_amount; 

    always_comb begin
        if (smaller) begin
            not_shifted = P_m_shifted;
            shifted = Z_m_padded;
        end else begin
            not_shifted = Z_m_padded;
            shifted = P_m_shifted;
        end
    end
endmodule
