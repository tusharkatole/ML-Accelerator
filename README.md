# ML-Accelerator(Non Maximum Suppression)

## CONTENTS
Final Codes
1. [A-IoU2(Our Design that eliminates multiplications and divisions)](#A-IoU2)
2. [A-IoU1(Standard IoU formulation, including multiplications and divisions)](#A-IoU1)
3. [Non Maximum Suppression using approximate IOU formula without pipelining](#Non-Maximum-Suppression-using-approximate-IOU-formula-without-pipelining)
4. [Non Maximum Suppression using general IOU formula without pipelining](#Non-Maximum-Suppression-using-general-IOU-formula-without-pipelining)
5. [Non Maximum Suppression using general IOU formula with pipelining](#Non-Maximum-Suppression-using-general-IOU-formula-with-pipelining)


# **A-IoU2(Our Design that eliminates multiplications and divisions)**

## A-IoU2

Our Design that eliminates multiplications and divisions

<details>
  <summary>Verilog: <code>calculate_iou</code> module (click to expand)</summary>


```
module calculate_iou (
    input wire clk,
    input wire reset,
    input wire start,
    input  [17:0] x1, y1, w1, h1,
    input  [17:0] x2, y2, w2, h2,
    output reg [17:0] iou,  // Approximate IoU in Q8.8 = (dx + dy) >> 8
    output reg done
);
    // Stage 0
    reg [17:0] x1_s0, y1_s0, w1_s0, h1_s0;
    reg [17:0] x2_s0, y2_s0, w2_s0, h2_s0;
    reg        start_s0;

    // Stage 1: Compute x+w, y+h
    reg [17:0] x1w1_s1, x2w2_s1, y1h1_s1, y2h2_s1;
    reg [17:0] x1_s1, x2_s1, y1_s1, y2_s1;
    reg        start_s1;

    // Stage 2: Compute dx, dy
    reg [17:0] dx_s2, dy_s2;
    reg        start_s2;

    // Stage 3: Output iou = (dx + dy) >> 8
    reg [17:0] sum_s3;
    reg        start_s3;

    always @(posedge clk) begin
        if (reset) begin
            iou <= 0;
            done <= 0;
            start_s0 <= 0;
            start_s1 <= 0;
            start_s2 <= 0;
            start_s3 <= 0;
        end else begin
            // Stage 0: Latch inputs
            if (start) begin
                x1_s0 <= x1; y1_s0 <= y1; w1_s0 <= w1; h1_s0 <= h1;
                x2_s0 <= x2; y2_s0 <= y2; w2_s0 <= w2; h2_s0 <= h2;
                start_s0 <= 1;
            end else begin
                start_s0 <= 0;
            end

            // Stage 1: Calculate box boundaries
            x1w1_s1 <= x1_s0 + w1_s0;
            x2w2_s1 <= x2_s0 + w2_s0;
            y1h1_s1 <= y1_s0 + h1_s0;
            y2h2_s1 <= y2_s0 + h2_s0;
            x1_s1 <= x1_s0; x2_s1 <= x2_s0;
            y1_s1 <= y1_s0; y2_s1 <= y2_s0;
            start_s1 <= start_s0;

            // Stage 2: Compute intersection dimensions
            dx_s2 <= (x1w1_s1 < x2w2_s1 ? x1w1_s1 : x2w2_s1) - (x1_s1 > x2_s1 ? x1_s1 : x2_s1);
            dy_s2 <= (y1h1_s1 < y2h2_s1 ? y1h1_s1 : y2h2_s1) - (y1_s1 > y2_s1 ? y1_s1 : y2_s1);
            start_s2 <= start_s1;

            // Stage 3: Final computation and output
            sum_s3 <= dx_s2 + dy_s2;
            start_s3 <= start_s2;

            // Output stage
            iou <= sum_s3 >> 8;  // convert to Q8.8
            done <= start_s3;
        end
    end

endmodule
```
</details>


<details>
  <summary>Verilog: <code>nms_top</code> module (click to expand)</summary>


## NMS_TOP_Code

```
module nms_top #(
    parameter NUM_BOXES = 10,                  // Number of boxes (configurable)
    parameter IOU_THRESHOLD = 18'd300,         // 0.7 in Q8.8 format
    parameter DATA_WIDTH = 18,
    parameter ADDR_WIDTH = $clog2(NUM_BOXES * 5), // Dynamic address width
    parameter TOTAL_ENTRIES = NUM_BOXES * 5    // Dynamic total entries
)(
    input wire clk,
    input wire reset,
    output wire [NUM_BOXES-1:0] suppressed
);

    // Block RAM signals
    wire [DATA_WIDTH-1:0] data_out;
    reg [ADDR_WIDTH-1:0] addr;
    reg [ADDR_WIDTH-1:0] delayed_addr_1, delayed_addr_2;
    
    // 2-cycle delay for BRAM output to track which address the data corresponds to
    always @(posedge clk) begin
        if (reset) begin
            delayed_addr_1 <= {ADDR_WIDTH{1'b1}}; // Invalid address initially
            delayed_addr_2 <= {ADDR_WIDTH{1'b1}}; // Invalid address initially
        end else begin
            delayed_addr_1 <= addr;
            delayed_addr_2 <= delayed_addr_1;
        end
    end

    // Block RAM instance
    blk_mem_gen_0 coe_rom (
        .clka(clk),
        .ena(1'b1),
        .wea(1'b0),
        .addra(addr),
        .dina({DATA_WIDTH{1'b0}}),
        .douta(data_out)
    );

    // Box storage registers - dynamically sized
    reg [DATA_WIDTH-1:0] boxes [0:NUM_BOXES-1][0:4]; // X, Y, W, H, S for each box
    reg [DATA_WIDTH-1:0] scores [0:NUM_BOXES-1];
    reg [NUM_BOXES-1:0] valid;
    reg [NUM_BOXES-1:0] suppressed_reg;
    
    assign suppressed = suppressed_reg;

    // FSM states
    localparam IDLE = 4'd0;
    localparam LOAD_BOXES = 4'd1;
    localparam WAIT_LAST_DATA = 4'd2;
    localparam SORT_INIT = 4'd3;
    localparam SORT_BOXES = 4'd4;
    localparam SORT_SWAP = 4'd5;
    localparam INIT_NMS = 4'd6;
    localparam COMPARE_BOXES = 4'd7;
    localparam CALCULATE_IOU = 4'd8;
    localparam WAIT_IOU = 4'd9;
    localparam UPDATE_VALID = 4'd10;
    localparam FINISH = 4'd11;
    
    reg [3:0] state;
    reg [ADDR_WIDTH-1:0] box_counter;
    reg [ADDR_WIDTH-1:0] current_box;
    reg [ADDR_WIDTH-1:0] compare_box;
    reg [2:0] wait_counter;
    reg [$clog2(NUM_BOXES):0] sort_passes; // Dynamic width for sort passes
    reg sorting_done;
    reg swap_needed;
    
    // IOU calculation signals
    wire [DATA_WIDTH-1:0] iou_result;
    reg iou_start;
    wire iou_done;
    reg [3:0] iou_wait_counter; // Counter for IOU pipeline delay - increased for 7-stage pipeline
    
    // Temporary registers for swapping
    reg [DATA_WIDTH-1:0] temp_score;
    reg [DATA_WIDTH-1:0] temp_box [0:4];
    
    // Registered inputs for IOU calculator for better timing
    reg [DATA_WIDTH-1:0] iou_x1, iou_y1, iou_w1, iou_h1;
    reg [DATA_WIDTH-1:0] iou_x2, iou_y2, iou_w2, iou_h2;
    
    // IOU calculator (optimized 7-stage pipeline)
    calculate_iou iou_calculator (
        .clk(clk),
        .reset(reset),
        .start(iou_start),
        .x1(iou_x1), .y1(iou_y1), .w1(iou_w1), .h1(iou_h1),
        .x2(iou_x2), .y2(iou_y2), .w2(iou_w2), .h2(iou_h2),
        .iou(iou_result),
        .done(iou_done)
    );

    // Synthesis-friendly initialization
    integer init_i, init_j;
    
    // Main FSM
    always @(posedge clk) begin
        if (reset) begin
            state <= IDLE;
            box_counter <= 0;
            current_box <= 0;
            compare_box <= 0;
            suppressed_reg <= {NUM_BOXES{1'b0}};
            valid <= {NUM_BOXES{1'b1}};
            addr <= 0;
            iou_start <= 0;
            wait_counter <= 0;
            sort_passes <= 0;
            sorting_done <= 0;
            swap_needed <= 0;
            temp_score <= 0;
            iou_wait_counter <= 0;
            
            // Initialize IOU input registers
            iou_x1 <= 0; iou_y1 <= 0; iou_w1 <= 0; iou_h1 <= 0;
            iou_x2 <= 0; iou_y2 <= 0; iou_w2 <= 0; iou_h2 <= 0;
            
            // Initialize arrays in synthesis-friendly way
            for (init_i = 0; init_i < NUM_BOXES; init_i = init_i + 1) begin
                scores[init_i] <= 0;
                for (init_j = 0; init_j < 5; init_j = init_j + 1) begin
                    boxes[init_i][init_j] <= 0;
                end
            end
            for (init_j = 0; init_j < 5; init_j = init_j + 1) begin
                temp_box[init_j] <= 0;
            end
        end else begin
            case (state)
                IDLE: begin
                    state <= LOAD_BOXES;
                    box_counter <= 0;
                    addr <= 0;
                    wait_counter <= 0;
                end
                
                LOAD_BOXES: begin
                    if (addr < TOTAL_ENTRIES) begin
                        addr <= addr + 1;
                    end else begin
                        // Wait for the last 2 data items to be processed
                        state <= WAIT_LAST_DATA;
                        wait_counter <= 0;
                    end
                    
                    // Store data using delayed_addr_2 (2-cycle delayed address)
                    // Only store when we have valid delayed address within range
                    if (delayed_addr_2 < TOTAL_ENTRIES && delayed_addr_2 != {ADDR_WIDTH{1'b1}}) begin
                        case (delayed_addr_2 % 5)
                            0: boxes[delayed_addr_2/5][0] <= data_out; // X
                            1: boxes[delayed_addr_2/5][1] <= data_out; // Y
                            2: boxes[delayed_addr_2/5][2] <= data_out; // W
                            3: boxes[delayed_addr_2/5][3] <= data_out; // H
                            4: begin
                                boxes[delayed_addr_2/5][4] <= data_out; // S
                                scores[delayed_addr_2/5] <= data_out;
                            end
                        endcase
                    end
                end
                
                WAIT_LAST_DATA: begin
                    // Wait 2 more cycles to get the last data from BRAM
                    if (wait_counter < 2) begin
                        wait_counter <= wait_counter + 1;
                        // Continue storing delayed data
                        if (delayed_addr_2 < TOTAL_ENTRIES && delayed_addr_2 != {ADDR_WIDTH{1'b1}}) begin
                            case (delayed_addr_2 % 5)
                                0: boxes[delayed_addr_2/5][0] <= data_out; // X
                                1: boxes[delayed_addr_2/5][1] <= data_out; // Y
                                2: boxes[delayed_addr_2/5][2] <= data_out; // W
                                3: boxes[delayed_addr_2/5][3] <= data_out; // H
                                4: begin
                                    boxes[delayed_addr_2/5][4] <= data_out; // S
                                    scores[delayed_addr_2/5] <= data_out;
                                end
                            endcase
                        end
                    end else begin
                        state <= SORT_INIT;
                        sort_passes <= 0;
                        sorting_done <= 0;
                    end
                end
                
                SORT_INIT: begin
                    box_counter <= 0;
                    state <= SORT_BOXES;
                end
                
                SORT_BOXES: begin
                    // Bubble sort implementation - one comparison per clock
                    if (box_counter < NUM_BOXES-1) begin
                        if (scores[box_counter] < scores[box_counter+1]) begin
                            // Store values in temp registers first
                            temp_score <= scores[box_counter];
                            temp_box[0] <= boxes[box_counter][0];
                            temp_box[1] <= boxes[box_counter][1];
                            temp_box[2] <= boxes[box_counter][2];
                            temp_box[3] <= boxes[box_counter][3];
                            temp_box[4] <= boxes[box_counter][4];
                            swap_needed <= 1'b1;
                            state <= SORT_SWAP;
                        end else begin
                            box_counter <= box_counter + 1;
                        end
                    end else begin
                        // One pass completed
                        sort_passes <= sort_passes + 1;
                        if (sort_passes < NUM_BOXES-1) begin
                            // Need more passes
                            state <= SORT_INIT;
                        end else begin
                            // Sorting complete
                            state <= INIT_NMS;
                            current_box <= 0;
                            suppressed_reg <= {NUM_BOXES{1'b0}};
                            valid <= {NUM_BOXES{1'b1}};
                        end
                    end
                end
                
                SORT_SWAP: begin
                    if (swap_needed) begin
                        // Perform the swap using temp registers
                        scores[box_counter] <= scores[box_counter+1];
                        scores[box_counter+1] <= temp_score;
                        boxes[box_counter][0] <= boxes[box_counter+1][0];
                        boxes[box_counter][1] <= boxes[box_counter+1][1];
                        boxes[box_counter][2] <= boxes[box_counter+1][2];
                        boxes[box_counter][3] <= boxes[box_counter+1][3];
                        boxes[box_counter][4] <= boxes[box_counter+1][4];
                        boxes[box_counter+1][0] <= temp_box[0];
                        boxes[box_counter+1][1] <= temp_box[1];
                        boxes[box_counter+1][2] <= temp_box[2];
                        boxes[box_counter+1][3] <= temp_box[3];
                        boxes[box_counter+1][4] <= temp_box[4];
                        swap_needed <= 1'b0;
                    end
                    box_counter <= box_counter + 1;
                    state <= SORT_BOXES;
                end
                
                INIT_NMS: begin
                    $display("DEBUG: INIT_NMS - current_box=%0d, NUM_BOXES=%0d, valid[%0d]=%b, suppressed[%0d]=%b", 
                             current_box, NUM_BOXES, current_box, 
                             (current_box < NUM_BOXES) ? valid[current_box] : 1'b0, 
                             current_box, 
                             (current_box < NUM_BOXES) ? suppressed_reg[current_box] : 1'b0);
                    if (current_box < NUM_BOXES) begin
                        if (valid[current_box] && !suppressed_reg[current_box]) begin
                            compare_box <= current_box + 1;
                            $display("DEBUG: Starting comparison for box %0d with boxes %0d onwards", current_box, current_box + 1);
                            state <= COMPARE_BOXES;
                        end else begin
                            $display("DEBUG: Skipping box %0d (valid=%b, suppressed=%b)", 
                                     current_box, valid[current_box], suppressed_reg[current_box]);
                            current_box <= current_box + 1;
                        end
                    end else begin
                        $display("DEBUG: NMS complete - all boxes processed");
                        state <= FINISH;
                    end
                end
                
                COMPARE_BOXES: begin
                    $display("DEBUG: COMPARE_BOXES - compare_box=%0d, NUM_BOXES=%0d, valid=%b, suppressed=%b", 
                             compare_box, NUM_BOXES, valid[compare_box], suppressed_reg[compare_box]);
                    if (compare_box < NUM_BOXES) begin
                        if (valid[compare_box] && !suppressed_reg[compare_box]) begin
                            $display("DEBUG: Calculating IOU between box %0d and box %0d", current_box, compare_box);
                            // Load IOU input registers
                            iou_x1 <= boxes[current_box][0];
                            iou_y1 <= boxes[current_box][1];
                            iou_w1 <= boxes[current_box][2];
                            iou_h1 <= boxes[current_box][3];
                            iou_x2 <= boxes[compare_box][0];
                            iou_y2 <= boxes[compare_box][1];
                            iou_w2 <= boxes[compare_box][2];
                            iou_h2 <= boxes[compare_box][3];
                            iou_start <= 1'b1;
                            state <= CALCULATE_IOU;
                        end else begin
                            $display("DEBUG: Skipping comparison with box %0d (valid=%b, suppressed=%b)", 
                                     compare_box, valid[compare_box], suppressed_reg[compare_box]);
                            compare_box <= compare_box + 1;
                            // Stay in COMPARE_BOXES state to check next box
                        end
                    end else begin
                        $display("DEBUG: Finished comparing box %0d with all others", current_box);
                        current_box <= current_box + 1;
                        state <= INIT_NMS;
                    end
                end
                
                CALCULATE_IOU: begin
                    iou_start <= 1'b0;
                    state <= WAIT_IOU;
                    iou_wait_counter <= 0;
                end
                
                WAIT_IOU: begin
                    // Wait for pipelined IOU calculation to complete
                    if (iou_done) begin
                        state <= UPDATE_VALID;
                    end else begin
                        // Additional safety counter to prevent infinite wait
                        if (iou_wait_counter < 10) begin // 9 pipeline stages + margin
                            iou_wait_counter <= iou_wait_counter + 1;
                        end else begin
                            state <= UPDATE_VALID; // Force progression
                        end
                    end
                end
                
                UPDATE_VALID: begin
                    $display("DEBUG: UPDATE_VALID - IOU=%0d, threshold=%0d, will suppress=%b", 
                             iou_result, IOU_THRESHOLD, (iou_result > IOU_THRESHOLD));
                    if (iou_result > IOU_THRESHOLD) begin
                        $display("DEBUG: Suppressing box %0d", compare_box);
                        suppressed_reg[compare_box] <= 1'b1;
                        valid[compare_box] <= 1'b0;
                    end else begin
                        $display("DEBUG: Keeping box %0d", compare_box);
                    end
                    compare_box <= compare_box + 1;
                    state <= COMPARE_BOXES;
                end
                
                FINISH: begin
                    // NMS complete
                    state <= FINISH;
                end
                
                default: begin
                    state <= IDLE;
                end
            endcase
        end
    end

   
 // Debug display for BRAM access
    always @(posedge clk) begin
        if (state == LOAD_BOXES || state == WAIT_LAST_DATA) begin
            if (delayed_addr_2 < TOTAL_ENTRIES && delayed_addr_2 != {ADDR_WIDTH{1'b1}}) begin
                $display("ADDR = %2d, Data = %6d", delayed_addr_2, data_out);
            end
        end
    end

    // Debug display for loaded boxes
    always @(posedge clk) begin
        if (state == SORT_INIT && sort_passes == 0) begin
            $display("=== Loaded Boxes (Before Sorting) ===");
            for (init_i = 0; init_i < NUM_BOXES; init_i = init_i + 1) begin
                $display("Box[%0d]: X=%0d, Y=%0d, W=%0d, H=%0d, S=%0d", 
                         init_i, boxes[init_i][0], boxes[init_i][1], boxes[init_i][2], boxes[init_i][3], boxes[init_i][4]);
            end
        end
    end

    // Debug display for sorted boxes
    always @(posedge clk) begin
        if (state == INIT_NMS && current_box == 0) begin
            $display("=== Sorted Boxes (By Score) ===");
            for (init_i = 0; init_i < NUM_BOXES; init_i = init_i + 1) begin
                $display("Box[%0d]: X=%0d, Y=%0d, W=%0d, H=%0d, S=%0d", 
                         init_i, boxes[init_i][0], boxes[init_i][1], boxes[init_i][2], boxes[init_i][3], boxes[init_i][4]);
            end
        end
    end

    // Debug display for NMS progress
    always @(posedge clk) begin
        if (state == INIT_NMS) begin
            $display("NMS: Processing current_box = %0d, valid = %b, suppressed = %b", 
                     current_box, valid, suppressed_reg);
        end
        if (state == COMPARE_BOXES) begin
            $display("NMS: Comparing current_box = %0d with compare_box = %0d", 
                     current_box, compare_box);
        end
    end

    // Debug display for IOU calculations
    always @(posedge clk) begin
        if (state == UPDATE_VALID) begin
            $display("IOU between Box[%0d] and Box[%0d] = %0d (threshold=%0d) -> %s", 
                     current_box, compare_box, iou_result, IOU_THRESHOLD,
                     (iou_result > IOU_THRESHOLD) ? "SUPPRESS" : "KEEP");
        end
    end

    // Debug display for final suppression results
    always @(posedge clk) begin
        if (state == FINISH) begin
            $display("=== NMS Results ===");
            $display("Suppressed mask: %b", suppressed_reg);
            for (init_i = 0; init_i < NUM_BOXES; init_i = init_i + 1) begin
                $display("Box[%0d]: %s (Score=%0d)", init_i, 
                         suppressed_reg[init_i] ? "SUPPRESSED" : "KEPT", scores[init_i]);
            end
            $display("NMS completed at time %0t", $time);
            $display("Suppression results: %b", suppressed_reg);
        end
    end

endmodule

```


</details>

# **A-IoU1(Standard IoU formulation, including multiplications and divisions)**

Standard IoU formulation, including multiplications and divisions

## A-IoU1

<details>
  <summary>Verilog: <code>calculate_iou</code> module (click to expand)</summary>

```
module calculate_iou (
    input wire clk,
    input wire reset,
    input wire start,
    input  [17:0] x1, y1, w1, h1,
    input  [17:0] x2, y2, w2, h2,
    output reg [17:0] iou,
    output reg done
);
    
    // Pipeline registers for stage 1
    reg [17:0] x1_reg, y1_reg, w1_reg, h1_reg;
    reg [17:0] x2_reg, y2_reg, w2_reg, h2_reg;
    reg [17:0] x1w1_reg, x2w2_reg, y1h1_reg, y2h2_reg;
    reg start_reg;
    
    // Pipeline registers for stage 2
    reg [17:0] dx_reg, dy_reg;
    reg [17:0] w1_reg2, h1_reg2, w2_reg2, h2_reg2;
    reg start_reg2;
    
    // Pipeline registers for stage 3
    reg [17:0] inter_reg;
    reg [17:0] denom_reg;
    reg start_reg3;
    
    // Combinational logic for each stage
    wire [17:0] x1w1_comb, x2w2_comb, y1h1_comb, y2h2_comb;
    wire [17:0] dx_comb, dy_comb;
    wire [17:0] inter_comb, denom_comb;
    wire [35:0] scaled_inter_comb;
    wire [17:0] iou_comb;
    
    // Stage 1: Calculate rectangle boundaries
    assign x1w1_comb = x1_reg + w1_reg;
    assign x2w2_comb = x2_reg + w2_reg;
    assign y1h1_comb = y1_reg + h1_reg;
    assign y2h2_comb = y2_reg + h2_reg;
    
    // Stage 2: Calculate intersection dimensions
    assign dx_comb = (x1w1_reg < x2w2_reg ? x1w1_reg : x2w2_reg) - (x1_reg > x2_reg ? x1_reg : x2_reg);
    assign dy_comb = (y1h1_reg < y2h2_reg ? y1h1_reg : y2h2_reg) - (y1_reg > y2_reg ? y1_reg : y2_reg);
    
    // Stage 3: Calculate intersection and union
    assign inter_comb = dx_reg + dy_reg;
    assign denom_comb = (w1_reg2 + h1_reg2 > w2_reg2 + h2_reg2) ? (w1_reg2 + h1_reg2) : (w2_reg2 + h2_reg2);
    
    // Stage 4: Scale and calculate final IOU
    assign scaled_inter_comb = (denom_reg != 0) ? (inter_reg * 256) : 0;
    assign iou_comb = (denom_reg != 0) ? scaled_inter_comb / denom_reg : 18'd0;
    
    // Pipeline control
    reg [2:0] pipeline_valid;
    
    always @(posedge clk) begin
        if (reset) begin
            // Reset all pipeline registers
            x1_reg <= 0; y1_reg <= 0; w1_reg <= 0; h1_reg <= 0;
            x2_reg <= 0; y2_reg <= 0; w2_reg <= 0; h2_reg <= 0;
            x1w1_reg <= 0; x2w2_reg <= 0; y1h1_reg <= 0; y2h2_reg <= 0;
            start_reg <= 0;
            
            dx_reg <= 0; dy_reg <= 0;
            w1_reg2 <= 0; h1_reg2 <= 0; w2_reg2 <= 0; h2_reg2 <= 0;
            start_reg2 <= 0;
            
            inter_reg <= 0; denom_reg <= 0;
            start_reg3 <= 0;
            
            iou <= 0;
            done <= 0;
            pipeline_valid <= 0;
        end else begin
            // Stage 1: Input registration and boundary calculation
            if (start) begin
                x1_reg <= x1; y1_reg <= y1; w1_reg <= w1; h1_reg <= h1;
                x2_reg <= x2; y2_reg <= y2; w2_reg <= w2; h2_reg <= h2;
                start_reg <= 1;
            end else begin
                start_reg <= 0;
            end
            
            // Stage 1 output registration
            x1w1_reg <= x1w1_comb;
            x2w2_reg <= x2w2_comb;
            y1h1_reg <= y1h1_comb;
            y2h2_reg <= y2h2_comb;
            
            // Stage 2: Intersection calculation
            dx_reg <= dx_comb;
            dy_reg <= dy_comb;
            w1_reg2 <= w1_reg;
            h1_reg2 <= h1_reg;
            w2_reg2 <= w2_reg;
            h2_reg2 <= h2_reg;
            start_reg2 <= start_reg;
            
            // Stage 3: Union calculation
            inter_reg <= inter_comb;
            denom_reg <= denom_comb;
            start_reg3 <= start_reg2;
            
            // Stage 4: Final IOU calculation
            iou <= iou_comb;
            done <= start_reg3;
            
            // Pipeline valid tracking
            pipeline_valid <= {pipeline_valid[1:0], start};
        end
    end

endmodule

```

</details>

<details>
  <summary>Verilog: <code>nms_top</code> module (click to expand)</summary>


## NMS_TOP_Code

```
module nms_top #(
    parameter NUM_BOXES = 10,                  // Number of boxes (configurable)
    parameter IOU_THRESHOLD = 18'd300,         // 0.7 in Q8.8 format
    parameter DATA_WIDTH = 18,
    parameter ADDR_WIDTH = $clog2(NUM_BOXES * 5), // Dynamic address width
    parameter TOTAL_ENTRIES = NUM_BOXES * 5    // Dynamic total entries
)(
    input wire clk,
    input wire reset,
    output wire [NUM_BOXES-1:0] suppressed
);

    // Block RAM signals
    wire [DATA_WIDTH-1:0] data_out;
    reg [ADDR_WIDTH-1:0] addr;
    reg [ADDR_WIDTH-1:0] delayed_addr_1, delayed_addr_2;
    
    // 2-cycle delay for BRAM output to track which address the data corresponds to
    always @(posedge clk) begin
        if (reset) begin
            delayed_addr_1 <= {ADDR_WIDTH{1'b1}}; // Invalid address initially
            delayed_addr_2 <= {ADDR_WIDTH{1'b1}}; // Invalid address initially
        end else begin
            delayed_addr_1 <= addr;
            delayed_addr_2 <= delayed_addr_1;
        end
    end

    // Block RAM instance
    blk_mem_gen_0 coe_rom (
        .clka(clk),
        .ena(1'b1),
        .wea(1'b0),
        .addra(addr),
        .dina({DATA_WIDTH{1'b0}}),
        .douta(data_out)
    );

    // Box storage registers - dynamically sized
    reg [DATA_WIDTH-1:0] boxes [0:NUM_BOXES-1][0:4]; // X, Y, W, H, S for each box
    reg [DATA_WIDTH-1:0] scores [0:NUM_BOXES-1];
    reg [NUM_BOXES-1:0] valid;
    reg [NUM_BOXES-1:0] suppressed_reg;
    
    assign suppressed = suppressed_reg;

    // FSM states
    localparam IDLE = 4'd0;
    localparam LOAD_BOXES = 4'd1;
    localparam WAIT_LAST_DATA = 4'd2;
    localparam SORT_INIT = 4'd3;
    localparam SORT_BOXES = 4'd4;
    localparam SORT_SWAP = 4'd5;
    localparam INIT_NMS = 4'd6;
    localparam COMPARE_BOXES = 4'd7;
    localparam CALCULATE_IOU = 4'd8;
    localparam WAIT_IOU = 4'd9;
    localparam UPDATE_VALID = 4'd10;
    localparam FINISH = 4'd11;
    
    reg [3:0] state;
    reg [ADDR_WIDTH-1:0] box_counter;
    reg [ADDR_WIDTH-1:0] current_box;
    reg [ADDR_WIDTH-1:0] compare_box;
    reg [2:0] wait_counter;
    reg [$clog2(NUM_BOXES):0] sort_passes; // Dynamic width for sort passes
    reg sorting_done;
    reg swap_needed;
    
    // IOU calculation signals
    wire [DATA_WIDTH-1:0] iou_result;
    reg iou_start;
    wire iou_done;
    reg [3:0] iou_wait_counter; // Counter for IOU pipeline delay - increased for 7-stage pipeline
    
    // Temporary registers for swapping
    reg [DATA_WIDTH-1:0] temp_score;
    reg [DATA_WIDTH-1:0] temp_box [0:4];
    
    // Registered inputs for IOU calculator for better timing
    reg [DATA_WIDTH-1:0] iou_x1, iou_y1, iou_w1, iou_h1;
    reg [DATA_WIDTH-1:0] iou_x2, iou_y2, iou_w2, iou_h2;
    
    // IOU calculator (optimized 7-stage pipeline)
    calculate_iou iou_calculator (
        .clk(clk),
        .reset(reset),
        .start(iou_start),
        .x1(iou_x1), .y1(iou_y1), .w1(iou_w1), .h1(iou_h1),
        .x2(iou_x2), .y2(iou_y2), .w2(iou_w2), .h2(iou_h2),
        .iou(iou_result),
        .done(iou_done)
    );

    // Synthesis-friendly initialization
    integer init_i, init_j;
    
    // Main FSM
    always @(posedge clk) begin
        if (reset) begin
            state <= IDLE;
            box_counter <= 0;
            current_box <= 0;
            compare_box <= 0;
            suppressed_reg <= {NUM_BOXES{1'b0}};
            valid <= {NUM_BOXES{1'b1}};
            addr <= 0;
            iou_start <= 0;
            wait_counter <= 0;
            sort_passes <= 0;
            sorting_done <= 0;
            swap_needed <= 0;
            temp_score <= 0;
            iou_wait_counter <= 0;
            
            // Initialize IOU input registers
            iou_x1 <= 0; iou_y1 <= 0; iou_w1 <= 0; iou_h1 <= 0;
            iou_x2 <= 0; iou_y2 <= 0; iou_w2 <= 0; iou_h2 <= 0;
            
            // Initialize arrays in synthesis-friendly way
            for (init_i = 0; init_i < NUM_BOXES; init_i = init_i + 1) begin
                scores[init_i] <= 0;
                for (init_j = 0; init_j < 5; init_j = init_j + 1) begin
                    boxes[init_i][init_j] <= 0;
                end
            end
            for (init_j = 0; init_j < 5; init_j = init_j + 1) begin
                temp_box[init_j] <= 0;
            end
        end else begin
            case (state)
                IDLE: begin
                    state <= LOAD_BOXES;
                    box_counter <= 0;
                    addr <= 0;
                    wait_counter <= 0;
                end
                
                LOAD_BOXES: begin
                    if (addr < TOTAL_ENTRIES) begin
                        addr <= addr + 1;
                    end else begin
                        // Wait for the last 2 data items to be processed
                        state <= WAIT_LAST_DATA;
                        wait_counter <= 0;
                    end
                    
                    // Store data using delayed_addr_2 (2-cycle delayed address)
                    // Only store when we have valid delayed address within range
                    if (delayed_addr_2 < TOTAL_ENTRIES && delayed_addr_2 != {ADDR_WIDTH{1'b1}}) begin
                        case (delayed_addr_2 % 5)
                            0: boxes[delayed_addr_2/5][0] <= data_out; // X
                            1: boxes[delayed_addr_2/5][1] <= data_out; // Y
                            2: boxes[delayed_addr_2/5][2] <= data_out; // W
                            3: boxes[delayed_addr_2/5][3] <= data_out; // H
                            4: begin
                                boxes[delayed_addr_2/5][4] <= data_out; // S
                                scores[delayed_addr_2/5] <= data_out;
                            end
                        endcase
                    end
                end
                
                WAIT_LAST_DATA: begin
                    // Wait 2 more cycles to get the last data from BRAM
                    if (wait_counter < 2) begin
                        wait_counter <= wait_counter + 1;
                        // Continue storing delayed data
                        if (delayed_addr_2 < TOTAL_ENTRIES && delayed_addr_2 != {ADDR_WIDTH{1'b1}}) begin
                            case (delayed_addr_2 % 5)
                                0: boxes[delayed_addr_2/5][0] <= data_out; // X
                                1: boxes[delayed_addr_2/5][1] <= data_out; // Y
                                2: boxes[delayed_addr_2/5][2] <= data_out; // W
                                3: boxes[delayed_addr_2/5][3] <= data_out; // H
                                4: begin
                                    boxes[delayed_addr_2/5][4] <= data_out; // S
                                    scores[delayed_addr_2/5] <= data_out;
                                end
                            endcase
                        end
                    end else begin
                        state <= SORT_INIT;
                        sort_passes <= 0;
                        sorting_done <= 0;
                    end
                end
                
                SORT_INIT: begin
                    box_counter <= 0;
                    state <= SORT_BOXES;
                end
                
                SORT_BOXES: begin
                    // Bubble sort implementation - one comparison per clock
                    if (box_counter < NUM_BOXES-1) begin
                        if (scores[box_counter] < scores[box_counter+1]) begin
                            // Store values in temp registers first
                            temp_score <= scores[box_counter];
                            temp_box[0] <= boxes[box_counter][0];
                            temp_box[1] <= boxes[box_counter][1];
                            temp_box[2] <= boxes[box_counter][2];
                            temp_box[3] <= boxes[box_counter][3];
                            temp_box[4] <= boxes[box_counter][4];
                            swap_needed <= 1'b1;
                            state <= SORT_SWAP;
                        end else begin
                            box_counter <= box_counter + 1;
                        end
                    end else begin
                        // One pass completed
                        sort_passes <= sort_passes + 1;
                        if (sort_passes < NUM_BOXES-1) begin
                            // Need more passes
                            state <= SORT_INIT;
                        end else begin
                            // Sorting complete
                            state <= INIT_NMS;
                            current_box <= 0;
                            suppressed_reg <= {NUM_BOXES{1'b0}};
                            valid <= {NUM_BOXES{1'b1}};
                        end
                    end
                end
                
                SORT_SWAP: begin
                    if (swap_needed) begin
                        // Perform the swap using temp registers
                        scores[box_counter] <= scores[box_counter+1];
                        scores[box_counter+1] <= temp_score;
                        boxes[box_counter][0] <= boxes[box_counter+1][0];
                        boxes[box_counter][1] <= boxes[box_counter+1][1];
                        boxes[box_counter][2] <= boxes[box_counter+1][2];
                        boxes[box_counter][3] <= boxes[box_counter+1][3];
                        boxes[box_counter][4] <= boxes[box_counter+1][4];
                        boxes[box_counter+1][0] <= temp_box[0];
                        boxes[box_counter+1][1] <= temp_box[1];
                        boxes[box_counter+1][2] <= temp_box[2];
                        boxes[box_counter+1][3] <= temp_box[3];
                        boxes[box_counter+1][4] <= temp_box[4];
                        swap_needed <= 1'b0;
                    end
                    box_counter <= box_counter + 1;
                    state <= SORT_BOXES;
                end
                
                INIT_NMS: begin
                    $display("DEBUG: INIT_NMS - current_box=%0d, NUM_BOXES=%0d, valid[%0d]=%b, suppressed[%0d]=%b", 
                             current_box, NUM_BOXES, current_box, 
                             (current_box < NUM_BOXES) ? valid[current_box] : 1'b0, 
                             current_box, 
                             (current_box < NUM_BOXES) ? suppressed_reg[current_box] : 1'b0);
                    if (current_box < NUM_BOXES) begin
                        if (valid[current_box] && !suppressed_reg[current_box]) begin
                            compare_box <= current_box + 1;
                            $display("DEBUG: Starting comparison for box %0d with boxes %0d onwards", current_box, current_box + 1);
                            state <= COMPARE_BOXES;
                        end else begin
                            $display("DEBUG: Skipping box %0d (valid=%b, suppressed=%b)", 
                                     current_box, valid[current_box], suppressed_reg[current_box]);
                            current_box <= current_box + 1;
                        end
                    end else begin
                        $display("DEBUG: NMS complete - all boxes processed");
                        state <= FINISH;
                    end
                end
                
                COMPARE_BOXES: begin
                    $display("DEBUG: COMPARE_BOXES - compare_box=%0d, NUM_BOXES=%0d, valid=%b, suppressed=%b", 
                             compare_box, NUM_BOXES, valid[compare_box], suppressed_reg[compare_box]);
                    if (compare_box < NUM_BOXES) begin
                        if (valid[compare_box] && !suppressed_reg[compare_box]) begin
                            $display("DEBUG: Calculating IOU between box %0d and box %0d", current_box, compare_box);
                            // Load IOU input registers
                            iou_x1 <= boxes[current_box][0];
                            iou_y1 <= boxes[current_box][1];
                            iou_w1 <= boxes[current_box][2];
                            iou_h1 <= boxes[current_box][3];
                            iou_x2 <= boxes[compare_box][0];
                            iou_y2 <= boxes[compare_box][1];
                            iou_w2 <= boxes[compare_box][2];
                            iou_h2 <= boxes[compare_box][3];
                            iou_start <= 1'b1;
                            state <= CALCULATE_IOU;
                        end else begin
                            $display("DEBUG: Skipping comparison with box %0d (valid=%b, suppressed=%b)", 
                                     compare_box, valid[compare_box], suppressed_reg[compare_box]);
                            compare_box <= compare_box + 1;
                            // Stay in COMPARE_BOXES state to check next box
                        end
                    end else begin
                        $display("DEBUG: Finished comparing box %0d with all others", current_box);
                        current_box <= current_box + 1;
                        state <= INIT_NMS;
                    end
                end
                
                CALCULATE_IOU: begin
                    iou_start <= 1'b0;
                    state <= WAIT_IOU;
                    iou_wait_counter <= 0;
                end
                
                WAIT_IOU: begin
                    // Wait for pipelined IOU calculation to complete
                    if (iou_done) begin
                        state <= UPDATE_VALID;
                    end else begin
                        // Additional safety counter to prevent infinite wait
                        if (iou_wait_counter < 10) begin // 9 pipeline stages + margin
                            iou_wait_counter <= iou_wait_counter + 1;
                        end else begin
                            state <= UPDATE_VALID; // Force progression
                        end
                    end
                end
                
                UPDATE_VALID: begin
                    $display("DEBUG: UPDATE_VALID - IOU=%0d, threshold=%0d, will suppress=%b", 
                             iou_result, IOU_THRESHOLD, (iou_result > IOU_THRESHOLD));
                    if (iou_result > IOU_THRESHOLD) begin
                        $display("DEBUG: Suppressing box %0d", compare_box);
                        suppressed_reg[compare_box] <= 1'b1;
                        valid[compare_box] <= 1'b0;
                    end else begin
                        $display("DEBUG: Keeping box %0d", compare_box);
                    end
                    compare_box <= compare_box + 1;
                    state <= COMPARE_BOXES;
                end
                
                FINISH: begin
                    // NMS complete
                    state <= FINISH;
                end
                
                default: begin
                    state <= IDLE;
                end
            endcase
        end
    end

   
 // Debug display for BRAM access
    always @(posedge clk) begin
        if (state == LOAD_BOXES || state == WAIT_LAST_DATA) begin
            if (delayed_addr_2 < TOTAL_ENTRIES && delayed_addr_2 != {ADDR_WIDTH{1'b1}}) begin
                $display("ADDR = %2d, Data = %6d", delayed_addr_2, data_out);
            end
        end
    end

    // Debug display for loaded boxes
    always @(posedge clk) begin
        if (state == SORT_INIT && sort_passes == 0) begin
            $display("=== Loaded Boxes (Before Sorting) ===");
            for (init_i = 0; init_i < NUM_BOXES; init_i = init_i + 1) begin
                $display("Box[%0d]: X=%0d, Y=%0d, W=%0d, H=%0d, S=%0d", 
                         init_i, boxes[init_i][0], boxes[init_i][1], boxes[init_i][2], boxes[init_i][3], boxes[init_i][4]);
            end
        end
    end

    // Debug display for sorted boxes
    always @(posedge clk) begin
        if (state == INIT_NMS && current_box == 0) begin
            $display("=== Sorted Boxes (By Score) ===");
            for (init_i = 0; init_i < NUM_BOXES; init_i = init_i + 1) begin
                $display("Box[%0d]: X=%0d, Y=%0d, W=%0d, H=%0d, S=%0d", 
                         init_i, boxes[init_i][0], boxes[init_i][1], boxes[init_i][2], boxes[init_i][3], boxes[init_i][4]);
            end
        end
    end

    // Debug display for NMS progress
    always @(posedge clk) begin
        if (state == INIT_NMS) begin
            $display("NMS: Processing current_box = %0d, valid = %b, suppressed = %b", 
                     current_box, valid, suppressed_reg);
        end
        if (state == COMPARE_BOXES) begin
            $display("NMS: Comparing current_box = %0d with compare_box = %0d", 
                     current_box, compare_box);
        end
    end

    // Debug display for IOU calculations
    always @(posedge clk) begin
        if (state == UPDATE_VALID) begin
            $display("IOU between Box[%0d] and Box[%0d] = %0d (threshold=%0d) -> %s", 
                     current_box, compare_box, iou_result, IOU_THRESHOLD,
                     (iou_result > IOU_THRESHOLD) ? "SUPPRESS" : "KEEP");
        end
    end

    // Debug display for final suppression results
    always @(posedge clk) begin
        if (state == FINISH) begin
            $display("=== NMS Results ===");
            $display("Suppressed mask: %b", suppressed_reg);
            for (init_i = 0; init_i < NUM_BOXES; init_i = init_i + 1) begin
                $display("Box[%0d]: %s (Score=%0d)", init_i, 
                         suppressed_reg[init_i] ? "SUPPRESSED" : "KEPT", scores[init_i]);
            end
            $display("NMS completed at time %0t", $time);
            $display("Suppression results: %b", suppressed_reg);
        end
    end

endmodule

```
</details>

# **Non Maximum Suppression using approximate IOU formula without pipelining**

## Pretreatment Unit (py)
<details>
  <summary>Python code to generate coordinates of Bounding Boxes</summary>

```c

import torch
import cv2
import numpy as np
from ultralytics import YOLO
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class YOLORawDetector:
    def __init__(self, model_path='yolov8n.pt'):
        """
        Initialize YOLO model for raw detection output
        Args:
            model_path: Path to YOLO model weights
        """
        self.model = YOLO(model_path)
        self.class_names = self.model.names
        
    def get_raw_detections(self, image_path, conf_threshold=0.01):
        """
        Get all raw detections without NMS
        Args:
            image_path: Path to input image
            conf_threshold: Minimum confidence threshold (very low to get all boxes)
        Returns:
            dict: Contains all detection information
        """
        # Load and process image
        image = cv2.imread(image_path)
        image_rgb = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
        
        # Run inference with very low confidence and no NMS
        results = self.model.predict(
            image_path,
            conf=conf_threshold,
            iou=1.0,  # Set IOU to 1.0 to disable NMS
            max_det=10000,  # Allow maximum detections
            verbose=False
        )
        
        # Extract raw detections
        detections = {
            'image_shape': image.shape,
            'total_detections': 0,
            'boxes': [],
            'confidences': [],
            'class_ids': [],
            'class_names': [],
            'all_detections': []
        }
        
        if len(results) > 0 and results[0].boxes is not None:
            boxes = results[0].boxes
            
            # Get all detection data
            if boxes.xyxy is not None:
                boxes_xyxy = boxes.xyxy.cpu().numpy()  # [x1, y1, x2, y2]
                confidences = boxes.conf.cpu().numpy()
                class_ids = boxes.cls.cpu().numpy().astype(int)
                
                detections['total_detections'] = len(boxes_xyxy)
                detections['boxes'] = boxes_xyxy
                detections['confidences'] = confidences
                detections['class_ids'] = class_ids
                detections['class_names'] = [self.class_names[cls_id] for cls_id in class_ids]
                
                # Create detailed detection list
                for i in range(len(boxes_xyxy)):
                    detection = {
                        'bbox': {
                            'x1': float(boxes_xyxy[i][0]),
                            'y1': float(boxes_xyxy[i][1]), 
                            'x2': float(boxes_xyxy[i][2]),
                            'y2': float(boxes_xyxy[i][3]),
                            'width': float(boxes_xyxy[i][2] - boxes_xyxy[i][0]),
                            'height': float(boxes_xyxy[i][3] - boxes_xyxy[i][1])
                        },
                        'confidence': float(confidences[i]),
                        'class_id': int(class_ids[i]),
                        'class_name': self.class_names[class_ids[i]]
                    }
                    detections['all_detections'].append(detection)
        
        return detections, image_rgb
    
    def print_all_detections(self, detections, min_conf=0.1):
        """
        Print all detections with details
        Args:
            detections: Detection dictionary from get_raw_detections
            min_conf: Minimum confidence to display
        """
        print(f"=== RAW YOLO DETECTIONS (No NMS) ===")
        print(f"Total detections found: {detections['total_detections']}")
        print(f"Image shape: {detections['image_shape']}")
        print("\n" + "="*80)
        
        # Sort by confidence (highest first)
        sorted_detections = sorted(detections['all_detections'], 
                                 key=lambda x: x['confidence'], reverse=True)
        
        displayed_count = 0
        for i, det in enumerate(sorted_detections):
            if det['confidence'] >= min_conf:
                print(f"Detection #{i+1}:")
                print(f"  Class: {det['class_name']} (ID: {det['class_id']})")
                print(f"  Confidence: {det['confidence']:.4f}")
                print(f"  Bounding Box:")
                print(f"    Bottom-left: ({det['bbox']['x1']:.1f}, {det['bbox']['y2']:.1f})")
                print(f"    Bottom-right: ({det['bbox']['x2']:.1f}, {det['bbox']['y2']:.1f})")
                print(f"    Width: {det['bbox']['width']:.1f}")
                print(f"    Height: {det['bbox']['height']:.1f}")
                print("-" * 50)
                displayed_count += 1
        
        print(f"\nDisplayed {displayed_count} detections with confidence >= {min_conf}")
        low_conf_count = len([d for d in sorted_detections if d['confidence'] < min_conf])
        if low_conf_count > 0:
            print(f"Hidden {low_conf_count} detections with confidence < {min_conf}")
    
    def visualize_detections(self, image, detections, min_conf=0.1, max_boxes=50):
        """
        Visualize detections on image
        Args:
            image: Input image (RGB format)
            detections: Detection dictionary
            min_conf: Minimum confidence to display
            max_boxes: Maximum number of boxes to display
        """
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        ax.imshow(image)
        
        # Filter and sort detections
        filtered_dets = [d for d in detections['all_detections'] if d['confidence'] >= min_conf]
        filtered_dets = sorted(filtered_dets, key=lambda x: x['confidence'], reverse=True)
        
        # Limit number of boxes for visualization
        display_dets = filtered_dets[:max_boxes]
        
        # Color map for different classes
        colors = plt.cm.tab10(np.linspace(0, 1, len(self.class_names)))
        
        for det in display_dets:
            bbox = det['bbox']
            class_id = det['class_id']
            conf = det['confidence']
            
            # Create rectangle
            rect = patches.Rectangle(
                (bbox['x1'], bbox['y1']), 
                bbox['width'], bbox['height'],
                linewidth=1, 
                edgecolor=colors[class_id % len(colors)], 
                facecolor='none',
                alpha=0.8
            )
            ax.add_patch(rect)
            
            # Add label
            label = f"{det['class_name']}: {conf:.3f}"
            ax.text(bbox['x1'], bbox['y1']-5, label, 
                   color=colors[class_id % len(colors)], 
                   fontsize=8, fontweight='bold',
                   bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.7))
        
        ax.set_title(f'Raw YOLO Detections (No NMS)\n'
                    f'Showing {len(display_dets)}/{len(filtered_dets)} detections '
                    f'(conf >= {min_conf})')
        ax.axis('off')
        plt.tight_layout()
        plt.show()
    
    def save_detections_to_file(self, detections, output_file='raw_detections.txt'):
        """
        Save all detections to text file
        """
        with open(output_file, 'w') as f:
            f.write("RAW YOLO DETECTIONS (No NMS)\n")
            f.write("="*50 + "\n")
            f.write(f"Total detections: {detections['total_detections']}\n")
            f.write(f"Image shape: {detections['image_shape']}\n\n")
            
            for i, det in enumerate(detections['all_detections']):
                f.write(f"Detection {i+1}:\n")
                f.write(f"  Class: {det['class_name']} (ID: {det['class_id']})\n")
                f.write(f"  Confidence: {det['confidence']:.6f}\n")
                f.write(f"  Bbox: [{det['bbox']['x1']:.2f}, {det['bbox']['y1']:.2f}, "
                       f"{det['bbox']['x2']:.2f}, {det['bbox']['y2']:.2f}]\n")
                f.write(f"  Size: {det['bbox']['width']:.2f} x {det['bbox']['height']:.2f}\n\n")
        
        print(f"Detections saved to {output_file}")

# Example usage
def main():
    # Initialize detector
    detector = YOLORawDetector('yolov8n.pt')  # or yolov8s.pt, yolov8m.pt, etc.
    
    # Replace with your image path
    image_path = '/kaggle/input/cat-dog/cat-dog1.jpg'
    
    try:
        # Get raw detections
        print("Processing image for raw detections...")
        detections, image = detector.get_raw_detections(
            image_path, 
            conf_threshold=0.01  # Very low threshold to get all boxes
        )
        
        # Print all detections
        detector.print_all_detections(detections, min_conf=0.1)
        
        # Visualize detections
        detector.visualize_detections(image, detections, min_conf=0.1, max_boxes=100)
        
        # Save to file
        detector.save_detections_to_file(detections)
        
        # Print summary statistics
        print(f"\n=== SUMMARY ===")
        print(f"Total raw detections: {detections['total_detections']}")
        
        # Confidence distribution
        confidences = detections['confidences']
        if len(confidences) > 0:
            print(f"Confidence range: {confidences.min():.4f} - {confidences.max():.4f}")
            print(f"Average confidence: {confidences.mean():.4f}")
            
            # Count by confidence levels
            high_conf = sum(c >= 0.5 for c in confidences)
            med_conf = sum(0.1 <= c < 0.5 for c in confidences)
            low_conf = sum(c < 0.1 for c in confidences)
            
            print(f"High confidence (≥0.5): {high_conf}")
            print(f"Medium confidence (0.1-0.5): {med_conf}")
            print(f"Low confidence (<0.1): {low_conf}")
        
    except Exception as e:
        print(f"Error: {e}")
        print("Make sure you have:")
        print("1. Installed ultralytics: pip install ultralytics")
        print("2. Provided a valid image path")
        print("3. Internet connection for first-time model download")

if __name__ == "__main__":
    main()

```

</details>


## NMS UNIT (py)
<details>
  <summary>Python for NMS using the coordinates generated from Pre Treatment Unit</summary>

```c
import cv2
import numpy as np
import matplotlib.pyplot as plt

# Your new detections (dog only)
detections = [
    ((248.3, 63.9, 380.1, 322.5), 0.9032, 16, "dog"),
    ((248.7, 64.8, 380.1, 323.1), 0.8971, 16, "dog"),
    ((249.0, 64.2, 380.6, 322.6), 0.8824, 16, "dog"),
    ((41.9, 62.9, 247.0, 323.0), 0.8742, 16, "dog"),
    ((246.8, 61.8, 380.0, 323.0), 0.8691, 16, "dog"),
    ((248.1, 64.4, 379.7, 322.4), 0.8653, 16, "dog"),
    ((248.8, 64.9, 381.3, 323.1), 0.8591, 16, "dog"),
    ((40.9, 62.3, 247.1, 323.0), 0.8586, 16, "dog"),
    ((41.3, 63.3, 248.1, 323.4), 0.8529, 16, "dog"),
    ((249.1, 64.9, 379.8, 322.7), 0.8252, 16, "dog"),
    ((40.5, 63.1, 247.3, 322.7), 0.8222, 16, "dog"),
    ((40.2, 61.9, 247.0, 323.0), 0.8093, 16, "dog"),
    ((250.3, 65.0, 381.8, 323.3), 0.7992, 16, "dog"),
    ((250.2, 64.9, 381.1, 323.5), 0.7875, 16, "dog"),
    ((39.4, 62.0, 247.9, 322.8), 0.7851, 16, "dog"),
    ((41.0, 63.0, 247.1, 323.0), 0.7465, 16, "dog"),
    ((41.0, 63.7, 248.5, 323.0), 0.7454, 16, "dog"),
    ((42.3, 63.1, 247.5, 323.4), 0.7267, 16, "dog"),
    ((251.2, 64.9, 378.6, 322.9), 0.6840, 16, "dog"),
    ((71.5, 63.9, 247.9, 324.9), 0.1707, 16, "dog"),
    ((39.9, 56.6, 248.4, 322.5), 0.1232, 16, "dog")
]

# NMS function
def nms(boxes, scores, threshold):
    if len(boxes) == 0:
        return [], []

    boxes = np.array(boxes)
    scores = np.array(scores)

    x1 = boxes[:, 0]
    y1 = boxes[:, 1]
    x2 = boxes[:, 2]
    y2 = boxes[:, 3]

    areas = (x2 - x1 + 1) * (y2 - y1 + 1)
    order = scores.argsort()[::-1]

    keep = []
    while order.size > 0:
        i = order[0]
        keep.append(i)

        xx1 = np.maximum(x1[i], x1[order[1:]])
        yy1 = np.maximum(y1[i], y1[order[1:]])
        xx2 = np.minimum(x2[i], x2[order[1:]])
        yy2 = np.minimum(y2[i], y2[order[1:]])

        w = np.maximum(0.0, xx2 - xx1 + 1)
        h = np.maximum(0.0, yy2 - yy1 + 1)
        inter = w * h

        iou = inter / (areas[i] + areas[order[1:]] - inter)
        inds = np.where(iou < threshold)[0]
        order = order[inds + 1]

    return [boxes[i] for i in keep], [scores[i] for i in keep]

# Prepare boxes and scores
boxes = [bbox for bbox, score, _, _ in detections]
scores = [score for _, score, _, _ in detections]
class_names = [class_name for _, _, _, class_name in detections]

# Apply NMS
nms_threshold = 0.5
final_boxes, final_scores = nms(boxes, scores, nms_threshold)

# Load image (replace with valid path)
image = cv2.imread("/kaggle/input/twodogs/twodog.jpg")
if image is None:
    raise FileNotFoundError("Image not found. Please check the path.")

# Make copies for visualization
original = image.copy()
nms_result = image.copy()

# Draw original boxes
for (x1, y1, x2, y2), score in zip(boxes, scores):
    cv2.rectangle(original, (int(x1), int(y1)), (int(x2), int(y2)), (0, 255, 255), 2)
    cv2.putText(original, f"dog:{score:.2f}", (int(x1), int(y1 - 5)),
                cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 0, 0), 1)

# Draw NMS-filtered boxes
for (x1, y1, x2, y2), score in zip(final_boxes, final_scores):
    cv2.rectangle(nms_result, (int(x1), int(y1)), (int(x2), int(y2)), (0, 255, 0), 2)
    cv2.putText(nms_result, f"dog:{score:.2f}", (int(x1), int(y1 - 5)),
                cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 0, 0), 1)

# Convert BGR to RGB for display
original_rgb = cv2.cvtColor(original, cv2.COLOR_BGR2RGB)
nms_rgb = cv2.cvtColor(nms_result, cv2.COLOR_BGR2RGB)

# Show both images
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.title("Original Detections")
plt.imshow(original_rgb)
plt.axis("off")

plt.subplot(1, 2, 2)
plt.title("After NMS (dog class)")
plt.imshow(nms_rgb)
plt.axis("off")

plt.tight_layout()
plt.show()

```
</details>


## IOU Calculation (py)
<details>
  <summary>Python code for iou calculation in fixed point specifically only for 2 given boxes</summary>

```c
def calculate_custom_iou_q8_8(x1, y1, w1, h1, x2, y2, w2, h2):
    x1w1 = x1 + w1
    x2w2 = x2 + w2
    y1h1 = y1 + h1
    y2h2 = y2 + h2

    dx = min(x1w1, x2w2) - max(x1, x2)
    dy = min(y1h1, y2h2) - max(y1, y2)

    dx = max(0, dx)
    dy = max(0, dy)

    inter = dx + dy
    denom = max(w1 + h1, w2 + h2)

    if denom != 0:
        scaled_inter = inter * 256
        iou_q8_8 = scaled_inter // denom  # Integer division
    else:
        iou_q8_8 = 0

    # Convert fixed Q8.8 to float for verification
    iou_float = iou_q8_8 / 256.0
    return iou_q8_8, iou_float

# Given values
x1, y1, w1, h1 = 7245, 83482, 81152, 42675
x2, y2, w2, h2 = 47206, 84582, 63462, 28314

iou_fixed, iou_float = calculate_custom_iou_q8_8(x1, y1, w1, h1, x2, y2, w2, h2)

print(f"IoU (Q8.8 fixed) : {iou_fixed} => binary: {bin(iou_fixed)}")
print(f"IoU (float)      : {iou_float:.6f}")

```
</details>


## NMS using floting point coordinates (py)
<details>
  <summary>Python code for NMS </summary>

```c
import numpy as np

# Input bounding boxes (x, y, h, w, s)
boxes = np.array([
    [28.3, 326.1, 317.0, 166.7, 0.8962],
    [27.5, 325.8, 316.0, 168.4, 0.8878],
    [28.5, 325.1, 314.5, 169.1, 0.8838],
    [27.1, 326.4, 314.9, 165.8, 0.834],
    [28.4, 324.1, 314.9, 169.7, 0.8328],
    [29.5, 324.8, 316.0, 166.8, 0.8238],
    [184.4, 330.4, 247.9, 110.6, 0.8174],
    [184.3, 330.2, 247.9, 110.5, 0.8128],
    [28.7, 324.0, 314.1, 172.9, 0.7984],
    [31.9, 325.0, 316.9, 172.7, 0.7921],
    [182.9, 331.0, 247.7, 112.0, 0.792],
    [185.3, 329.1, 245.9, 111.0, 0.7497],
    [32.6, 325.5, 317.6, 171.4, 0.748],
    [28.9, 327.0, 316.5, 164.1, 0.7129],
    [183.6, 329.2, 245.7, 112.4, 0.682],
    [179.3, 331.7, 246.9, 114.2, 0.6441],
    [188.4, 327.1, 243.5, 102.4, 0.5816],
    [176.7, 331.9, 247.1, 116.9, 0.5293],
    [25.4, 329.2, 272.2, 270.5, 0.28],
    [134.9, 328.6, 264.1, 160.3, 0.2092],
    [49.0, 328.7, 269.6, 246.9, 0.2069],
    [20.9, 329.2, 270.1, 275.0, 0.1998],
    [187.8, 327.2, 243.6, 104.5, 0.1751],
    [23.9, 328.9, 283.1, 273.3, 0.1],
])

# Modified IoU function
def modified_iou(box1, box2):
    x1, y1, h1, w1 = box1[:4]
    x2, y2, h2, w2 = box2[:4]

    x_overlap = max(0, min(x1 + w1, x2 + w2) - max(x1, x2))
    y_overlap = max(0, min(y1 + h1, y2 + h2) - max(y1, y2))
    numerator = x_overlap + y_overlap
    denominator = max(h1 + w1, h2 + w2)

    return numerator / denominator if denominator != 0 else 0

# Non-Maximum Suppression using modified IoU
def nms_modified(boxes, iou_threshold=0.5):
    # Sort boxes by score descending
    boxes = boxes[boxes[:, 4].argsort()[::-1]]
    keep = []

    while len(boxes) > 0:
        current = boxes[0]
        keep.append(current)
        remaining = []

        for box in boxes[1:]:
            iou = modified_iou(current, box)
            if iou < iou_threshold:
                remaining.append(box)

        boxes = np.array(remaining)

    return np.array(keep)

# Run NMS
result_boxes = nms_modified(boxes, iou_threshold=0.7)

# Print remaining unsuppressed boxes
print("Remaining bounding boxes (x, y, h, w, s):")
for b in result_boxes:
    print(f"{b[0]:.1f}, {b[1]:.1f}, {b[2]:.1f}, {b[3]:.1f}, {b[4]:.4f}")

```
</details>


## IOU Calculation Modified Algorithm (Verilog)
<details>
  <summary>This is the sub moudle for NMS UNIT</summary>

```c
module calculate_iou (
    //input wire clk,
    //input wire reset,
    //input wire start,
    input  [17:0] x1, y1, w1, h1,
    input  [17:0] x2, y2, w2, h2,
    output [17:0] iou
    //output reg done
);
    wire [17:0] x1w1 = x1 + w1;
    wire [17:0] x2w2 = x2 + w2;
    wire [17:0] y1h1 = y1 + h1;
    wire [17:0] y2h2 = y2 + h2;

    wire [17:0] dx = (x1w1 < x2w2 ? x1w1 : x2w2) - (x1 > x2 ? x1 : x2);
    wire [17:0] dy = (y1h1 < y2h2 ? y1h1 : y2h2) - (y1 > y2 ? y1 : y2);

    wire [17:0] inter = dx + dy;
    wire [17:0] denom = (w1 + h1 > w2 + h2) ? (w1 + h1) : (w2 + h2);

 
    wire [35:0] scaled_inter = (denom != 0) ? (inter * 256) : 0;
    assign iou = (denom != 0) ? scaled_inter / denom : 18'd0;

endmodule
```
</details>


## NMS UNIT (Verilog)
<details>
  <summary>Compelete top module for NMS</summary>

```c
module nms_top #(
    parameter NUM_BOXES = 8,                  // 8 boxes
    parameter IOU_THRESHOLD = 18'd180,        // 0.7 in Q8.8 format
    parameter DATA_WIDTH = 18,
    parameter ADDR_WIDTH = 6,                 // 6 bits for 8 boxes (8*5=40 needs 6 bits)
    parameter TOTAL_ENTRIES = 40              // 8 boxes * 5 parameters each
)(
    input wire clk,
    input wire reset,
    output wire [NUM_BOXES-1:0] suppressed
);

    // Block RAM signals
    wire [DATA_WIDTH-1:0] data_out;
    reg [ADDR_WIDTH-1:0] addr;
    reg [ADDR_WIDTH-1:0] delayed_addr_1, delayed_addr_2;
    
    // 2-cycle delay for BRAM output to track which address the data corresponds to
    always @(posedge clk) begin
        if (reset) begin
            delayed_addr_1 <= 6'b111111; // Invalid address initially
            delayed_addr_2 <= 6'b111111; // Invalid address initially
        end else begin
            delayed_addr_1 <= addr;
            delayed_addr_2 <= delayed_addr_1;
        end
    end

    // Block RAM instance
    blk_mem_gen_0 coe_rom (
        .clka(clk),
        .ena(1'b1),
        .wea(1'b0),
        .addra(addr),
        .dina({DATA_WIDTH{1'b0}}),
        .douta(data_out)
    );

    // Box storage registers - sized for 8 boxes
    reg [DATA_WIDTH-1:0] boxes [0:NUM_BOXES-1][0:4]; // X, Y, W, H, S for each box
    reg [DATA_WIDTH-1:0] scores [0:NUM_BOXES-1];
    reg [NUM_BOXES-1:0] valid;
    reg [NUM_BOXES-1:0] suppressed_reg;
    
    assign suppressed = suppressed_reg;

    // FSM states
    localparam IDLE = 4'd0;
    localparam LOAD_BOXES = 4'd1;
    localparam WAIT_LAST_DATA = 4'd2;
    localparam SORT_INIT = 4'd3;
    localparam SORT_BOXES = 4'd4;
    localparam SORT_SWAP = 4'd5;
    localparam INIT_NMS = 4'd6;
    localparam COMPARE_BOXES = 4'd7;
    localparam CALCULATE_IOU = 4'd8;
    localparam UPDATE_VALID = 4'd9;
    localparam FINISH = 4'd10;
    
    reg [3:0] state;
    reg [ADDR_WIDTH-1:0] box_counter;
    reg [ADDR_WIDTH-1:0] current_box;
    reg [ADDR_WIDTH-1:0] compare_box;
    reg [2:0] wait_counter;
    reg [3:0] sort_passes;
    reg sorting_done;
    reg swap_needed;
    
    // IOU calculation signals
    wire [DATA_WIDTH-1:0] iou_result;
    reg iou_start;
    wire iou_done = 1'b1; // Combinational - always done
    integer i, j;
    
    // Temporary registers for swapping
    reg [DATA_WIDTH-1:0] temp_score;
    reg [DATA_WIDTH-1:0] temp_box [0:4];
    
    // IOU calculator (combinational)
    calculate_iou iou_calculator (
        .x1(boxes[current_box][0]), .y1(boxes[current_box][1]), 
        .w1(boxes[current_box][2]), .h1(boxes[current_box][3]),
        .x2(boxes[compare_box][0]), .y2(boxes[compare_box][1]), 
        .w2(boxes[compare_box][2]), .h2(boxes[compare_box][3]),
        .iou(iou_result)
    );

    // Main FSM
    always @(posedge clk or posedge reset) begin
        if (reset) begin
            state <= IDLE;
            box_counter <= 0;
            current_box <= 0;
            compare_box <= 0;
            suppressed_reg <= {NUM_BOXES{1'b0}};
            valid <= {NUM_BOXES{1'b1}};
            addr <= 0;
            iou_start <= 0;
            wait_counter <= 0;
            sort_passes <= 0;
            sorting_done <= 0;
            swap_needed <= 0;
            
            // Initialize arrays
            for (i = 0; i < NUM_BOXES; i = i + 1) begin
                scores[i] <= 0;
                for (j = 0; j < 5; j = j + 1) begin
                    boxes[i][j] <= 0;
                end
                for (j = 0; j < 5; j = j + 1) begin
                    temp_box[j] <= 0;
                end
            end
            temp_score <= 0;
        end else begin
            case (state)
                IDLE: begin
                    state <= LOAD_BOXES;
                    box_counter <= 0;
                    addr <= 0;
                    wait_counter <= 0;
                end
                
                LOAD_BOXES: begin
                    if (addr < TOTAL_ENTRIES) begin  // 0 to 39 for 8 boxes
                        addr <= addr + 1;
                    end else begin
                        // Wait for the last 2 data items to be processed
                        state <= WAIT_LAST_DATA;
                        wait_counter <= 0;
                    end
                    
                    // Store data using delayed_addr_2 (2-cycle delayed address)
                    // Only store when we have valid delayed address within range
                    if (delayed_addr_2 < TOTAL_ENTRIES && delayed_addr_2 != 6'b111111) begin
                        case (delayed_addr_2 % 5)
                            0: boxes[delayed_addr_2/5][0] <= data_out; // X
                            1: boxes[delayed_addr_2/5][1] <= data_out; // Y
                            2: boxes[delayed_addr_2/5][2] <= data_out; // W
                            3: boxes[delayed_addr_2/5][3] <= data_out; // H
                            4: begin
                                boxes[delayed_addr_2/5][4] <= data_out; // S
                                scores[delayed_addr_2/5] <= data_out;
                            end
                        endcase
                    end
                end
                
                WAIT_LAST_DATA: begin
                    // Wait 2 more cycles to get the last data from BRAM
                    if (wait_counter < 2) begin
                        wait_counter <= wait_counter + 1;
                        // Continue storing delayed data
                        if (delayed_addr_2 < TOTAL_ENTRIES && delayed_addr_2 != 6'b111111) begin
                            case (delayed_addr_2 % 5)
                                0: boxes[delayed_addr_2/5][0] <= data_out; // X
                                1: boxes[delayed_addr_2/5][1] <= data_out; // Y
                                2: boxes[delayed_addr_2/5][2] <= data_out; // W
                                3: boxes[delayed_addr_2/5][3] <= data_out; // H
                                4: begin
                                    boxes[delayed_addr_2/5][4] <= data_out; // S
                                    scores[delayed_addr_2/5] <= data_out;
                                end
                            endcase
                        end
                    end else begin
                        state <= SORT_INIT;
                        sort_passes <= 0;
                        sorting_done <= 0;
                    end
                end
                
                SORT_INIT: begin
                    box_counter <= 0;
                    state <= SORT_BOXES;
                end
                
                SORT_BOXES: begin
                    // Bubble sort implementation - one comparison per clock
                    if (box_counter < NUM_BOXES-1) begin
                        if (scores[box_counter] < scores[box_counter+1]) begin
                            // Store values in temp registers first
                            temp_score <= scores[box_counter];
                            for (i = 0; i < 5; i = i + 1) begin
                                temp_box[i] <= boxes[box_counter][i];
                            end
                            swap_needed <= 1'b1;
                            state <= SORT_SWAP;
                        end else begin
                            box_counter <= box_counter + 1;
                        end
                    end else begin
                        // One pass completed
                        sort_passes <= sort_passes + 1;
                        if (sort_passes < NUM_BOXES-1) begin
                            // Need more passes
                            state <= SORT_INIT;
                        end else begin
                            // Sorting complete
                            state <= INIT_NMS;
                            current_box <= 0;
                            suppressed_reg <= {NUM_BOXES{1'b0}};
                            valid <= {NUM_BOXES{1'b1}};
                        end
                    end
                end
                
                SORT_SWAP: begin
                    if (swap_needed) begin
                        // Perform the swap using temp registers
                        scores[box_counter] <= scores[box_counter+1];
                        scores[box_counter+1] <= temp_score;
                        for (i = 0; i < 5; i = i + 1) begin
                            boxes[box_counter][i] <= boxes[box_counter+1][i];
                            boxes[box_counter+1][i] <= temp_box[i];
                        end
                        swap_needed <= 1'b0;
                    end
                    box_counter <= box_counter + 1;
                    state <= SORT_BOXES;
                end
                
                INIT_NMS: begin
                    $display("DEBUG: INIT_NMS - current_box=%0d, NUM_BOXES=%0d, valid[%0d]=%b, suppressed[%0d]=%b", 
                             current_box, NUM_BOXES, current_box, 
                             (current_box < NUM_BOXES) ? valid[current_box] : 1'b0, 
                             current_box, 
                             (current_box < NUM_BOXES) ? suppressed_reg[current_box] : 1'b0);
                    if (current_box < NUM_BOXES) begin
                        if (valid[current_box] && !suppressed_reg[current_box]) begin
                            compare_box <= current_box + 1;
                            $display("DEBUG: Starting comparison for box %0d with boxes %0d onwards", current_box, current_box + 1);
                            state <= COMPARE_BOXES;
                        end else begin
                            $display("DEBUG: Skipping box %0d (valid=%b, suppressed=%b)", 
                                     current_box, valid[current_box], suppressed_reg[current_box]);
                            current_box <= current_box + 1;
                        end
                    end else begin
                        $display("DEBUG: NMS complete - all boxes processed");
                        state <= FINISH;
                    end
                end
                
                COMPARE_BOXES: begin
                    $display("DEBUG: COMPARE_BOXES - compare_box=%0d, NUM_BOXES=%0d, valid=%b, suppressed=%b", 
                             compare_box, NUM_BOXES, valid[compare_box], suppressed_reg[compare_box]);
                    if (compare_box < NUM_BOXES) begin
                        if (valid[compare_box] && !suppressed_reg[compare_box]) begin
                            $display("DEBUG: Calculating IOU between box %0d and box %0d", current_box, compare_box);
                            iou_start <= 1'b1;
                            state <= CALCULATE_IOU;
                        end else begin
                            $display("DEBUG: Skipping comparison with box %0d (valid=%b, suppressed=%b)", 
                                     compare_box, valid[compare_box], suppressed_reg[compare_box]);
                            compare_box <= compare_box + 1;
                            // Stay in COMPARE_BOXES state to check next box
                        end
                    end else begin
                        $display("DEBUG: Finished comparing box %0d with all others", current_box);
                        current_box <= current_box + 1;
                        state <= INIT_NMS;
                    end
                end
                
                CALCULATE_IOU: begin
                    iou_start <= 1'b0;
                    state <= UPDATE_VALID; // Combinational IOU - go directly to update
                end
                
                UPDATE_VALID: begin
                    $display("DEBUG: UPDATE_VALID - IOU=%0d, threshold=%0d, will suppress=%b", 
                             iou_result, IOU_THRESHOLD, (iou_result > IOU_THRESHOLD));
                    if (iou_result > IOU_THRESHOLD) begin
                        $display("DEBUG: Suppressing box %0d", compare_box);
                        suppressed_reg[compare_box] <= 1'b1;
                        valid[compare_box] <= 1'b0;
                    end else begin
                        $display("DEBUG: Keeping box %0d", compare_box);
                    end
                    compare_box <= compare_box + 1;
                    state <= COMPARE_BOXES;
                end
                
                FINISH: begin
                    // NMS complete
                    state <= FINISH;
                end
                
                default: begin
                    state <= IDLE;
                end
            endcase
        end
    end

    // Debug display for BRAM access
    always @(posedge clk) begin
        if (state == LOAD_BOXES || state == WAIT_LAST_DATA) begin
            if (delayed_addr_2 < TOTAL_ENTRIES && delayed_addr_2 != 6'b111111) begin
                $display("ADDR = %2d, Data = %6d", delayed_addr_2, data_out);
            end
        end
    end

    // Debug display for loaded boxes
    always @(posedge clk) begin
        if (state == SORT_INIT && sort_passes == 0) begin
            $display("=== Loaded Boxes (Before Sorting) ===");
            for (i = 0; i < NUM_BOXES; i = i + 1) begin
                $display("Box[%0d]: X=%0d, Y=%0d, W=%0d, H=%0d, S=%0d", 
                         i, boxes[i][0], boxes[i][1], boxes[i][2], boxes[i][3], boxes[i][4]);
            end
        end
    end

    // Debug display for sorted boxes
    always @(posedge clk) begin
        if (state == INIT_NMS && current_box == 0) begin
            $display("=== Sorted Boxes (By Score) ===");
            for (i = 0; i < NUM_BOXES; i = i + 1) begin
                $display("Box[%0d]: X=%0d, Y=%0d, W=%0d, H=%0d, S=%0d", 
                         i, boxes[i][0], boxes[i][1], boxes[i][2], boxes[i][3], boxes[i][4]);
            end
        end
    end

    // Debug display for NMS progress
    always @(posedge clk) begin
        if (state == INIT_NMS) begin
            $display("NMS: Processing current_box = %0d, valid = %b, suppressed = %b", 
                     current_box, valid, suppressed_reg);
        end
        if (state == COMPARE_BOXES) begin
            $display("NMS: Comparing current_box = %0d with compare_box = %0d", 
                     current_box, compare_box);
        end
    end

    // Debug display for IOU calculations
    always @(posedge clk) begin
        if (state == UPDATE_VALID) begin
            $display("IOU between Box[%0d] and Box[%0d] = %0d (threshold=%0d) -> %s", 
                     current_box, compare_box, iou_result, IOU_THRESHOLD,
                     (iou_result > IOU_THRESHOLD) ? "SUPPRESS" : "KEEP");
        end
    end

    // Debug display for final suppression results
    always @(posedge clk) begin
        if (state == FINISH) begin
            $display("=== NMS Results ===");
            $display("Suppressed mask: %b", suppressed_reg);
            for (i = 0; i < NUM_BOXES; i = i + 1) begin
                $display("Box[%0d]: %s (Score=%0d)", i, 
                         suppressed_reg[i] ? "SUPPRESSED" : "KEPT", scores[i]);
            end
            $display("NMS completed at time %0t", $time);
            $display("Suppression results: %b", suppressed_reg);
        end
    end

endmodule

```
</details>

## Testbench (Verilog)
<details>
  <summary>All necessary signals are provided in testbench</summary>

```c
`timescale 1ns / 1ps

module nms_top_tb();

    // Parameters
    parameter NUM_BOXES = 8;
    parameter CLOCK_PERIOD = 10; // 10ns = 100MHz

    // Inputs
    reg clk;
    reg reset;

    // Outputs
    wire [NUM_BOXES-1:0] suppressed;

    // Instantiate the Unit Under Test (UUT)
    nms_top #(
        .NUM_BOXES(NUM_BOXES),
        .IOU_THRESHOLD(18'd180)  // 0.7 in Q8.8 format
    ) uut (
        .clk(clk),
        .reset(reset),
        .suppressed(suppressed)
    );

    // Clock generation
    initial begin
        clk = 0;
        forever #(CLOCK_PERIOD/2) clk = ~clk;
    end

    // Reset and stimulus
    initial begin
        // Initialize and reset
        reset = 1;
        #100;
        reset = 0;
        
        $display("Starting NMS simulation...");
        $display("Design will automatically read from block RAM");
        
        // Wait for completion (state == FINISH)
        wait(uut.state == 3'd7);
        #1000;
        
        // Display results
//        $display("NMS completed at time %t ns", $time);
//        $display("Suppression results: %b", suppressed);
        
//        // Expected output check
//        if (suppressed === 5'b01110) begin
//            $display("TEST PASSED: Correct suppression pattern detected");
//        end else begin
//            $display("TEST FAILED: Expected 01110, got %b", suppressed);
//        end
        
        
        $finish;
    end
endmodule
```
</details>



# Non Maximum Suppression using general IOU formula without pipelining

## IOU Calculation general Algorithm (Verilog)

<details>
  <summary>This is the sub moudle for NMS UNIT</summary>

```c
module calculate_iou (
    input  [17:0] x1, y1, w1, h1,
    input  [17:0] x2, y2, w2, h2,
    output [17:0] iou
);
    // Calculate box boundaries
    wire [17:0] x1w1 = x1 + w1;
    wire [17:0] x2w2 = x2 + w2;
    wire [17:0] y1h1 = y1 + h1;
    wire [17:0] y2h2 = y2 + h2;

    // Intersection coordinates
    wire [17:0] inter_x1 = (x1 > x2) ? x1 : x2;
    wire [17:0] inter_y1 = (y1 > y2) ? y1 : y2;
    wire [17:0] inter_x2 = (x1w1 < x2w2) ? x1w1 : x2w2;
    wire [17:0] inter_y2 = (y1h1 < y2h2) ? y1h1 : y2h2;

    // Intersection dimensions (clamped to prevent negatives)
    wire [17:0] inter_w = (inter_x2 > inter_x1) ? (inter_x2 - inter_x1) : 0;
    wire [17:0] inter_h = (inter_y2 > inter_y1) ? (inter_y2 - inter_y1) : 0;

    // Area calculations
    wire [35:0] area_inter = inter_w * inter_h;
    wire [35:0] area1 = w1 * h1;
    wire [35:0] area2 = w2 * h2;
    wire [35:0] area_union = area1 + area2 - area_inter;

    // IoU calculation with fixed-point scaling
    wire [43:0] scaled_iou = (area_union != 0) ? (area_inter * 256) / area_union : 0;
    assign iou = scaled_iou[17:0];
endmodule

```
</details>


## NMS UNIT (Verilog)
<details>
  <summary>Compelete top module for NMS</summary>

```c
module nms_top #(
    parameter NUM_BOXES = 11,                  // Number of boxes (configurable)
    parameter IOU_THRESHOLD = 18'd180,         // 0.7 in Q8.8 format
    parameter DATA_WIDTH = 18,
    parameter ADDR_WIDTH = $clog2(NUM_BOXES * 5), // Dynamic address width
    parameter TOTAL_ENTRIES = NUM_BOXES * 5    // Dynamic total entries
)(
    input wire clk,
    input wire reset,
    output wire [NUM_BOXES-1:0] suppressed
);

    // Block RAM signals
    wire [DATA_WIDTH-1:0] data_out;
    reg [ADDR_WIDTH-1:0] addr;
    reg [ADDR_WIDTH-1:0] delayed_addr_1, delayed_addr_2;
    
    // 2-cycle delay for BRAM output to track which address the data corresponds to
    always @(posedge clk) begin
        if (reset) begin
            delayed_addr_1 <= {ADDR_WIDTH{1'b1}}; // Invalid address initially
            delayed_addr_2 <= {ADDR_WIDTH{1'b1}}; // Invalid address initially
        end else begin
            delayed_addr_1 <= addr;
            delayed_addr_2 <= delayed_addr_1;
        end
    end

    // Block RAM instance
    blk_mem_gen_0 coe_rom (
        .clka(clk),
        .ena(1'b1),
        .wea(1'b0),
        .addra(addr),
        .dina({DATA_WIDTH{1'b0}}),
        .douta(data_out)
    );

    // Box storage registers - dynamically sized
    reg [DATA_WIDTH-1:0] boxes [0:NUM_BOXES-1][0:4]; // X, Y, W, H, S for each box
    reg [DATA_WIDTH-1:0] scores [0:NUM_BOXES-1];
    reg [NUM_BOXES-1:0] valid;
    reg [NUM_BOXES-1:0] suppressed_reg;
    
    assign suppressed = suppressed_reg;

    // FSM states
    localparam IDLE = 4'd0;
    localparam LOAD_BOXES = 4'd1;
    localparam WAIT_LAST_DATA = 4'd2;
    localparam SORT_INIT = 4'd3;
    localparam SORT_BOXES = 4'd4;
    localparam SORT_SWAP = 4'd5;
    localparam INIT_NMS = 4'd6;
    localparam COMPARE_BOXES = 4'd7;
    localparam CALCULATE_IOU = 4'd8;
    localparam UPDATE_VALID = 4'd9;
    localparam FINISH = 4'd10;
    
    reg [3:0] state;
    reg [ADDR_WIDTH-1:0] box_counter;
    reg [ADDR_WIDTH-1:0] current_box;
    reg [ADDR_WIDTH-1:0] compare_box;
    reg [2:0] wait_counter;
    reg [$clog2(NUM_BOXES):0] sort_passes; // Dynamic width for sort passes
    reg sorting_done;
    reg swap_needed;
    
    // IOU calculation signals
    wire [DATA_WIDTH-1:0] iou_result;
    reg iou_start;
    wire iou_done = 1'b1; // Combinational - always done
    
    // Temporary registers for swapping
    reg [DATA_WIDTH-1:0] temp_score;
    reg [DATA_WIDTH-1:0] temp_box [0:4];
    
    // IOU calculator (combinational)
    calculate_iou iou_calculator (
        .x1(boxes[current_box][0]), .y1(boxes[current_box][1]), 
        .w1(boxes[current_box][2]), .h1(boxes[current_box][3]),
        .x2(boxes[compare_box][0]), .y2(boxes[compare_box][1]), 
        .w2(boxes[compare_box][2]), .h2(boxes[compare_box][3]),
        .iou(iou_result)
    );

    // Synthesis-friendly initialization
    integer init_i, init_j;
    
    // Main FSM
    always @(posedge clk) begin
        if (reset) begin
            state <= IDLE;
            box_counter <= 0;
            current_box <= 0;
            compare_box <= 0;
            suppressed_reg <= {NUM_BOXES{1'b0}};
            valid <= {NUM_BOXES{1'b1}};
            addr <= 0;
            iou_start <= 0;
            wait_counter <= 0;
            sort_passes <= 0;
            sorting_done <= 0;
            swap_needed <= 0;
            temp_score <= 0;
            
            // Initialize arrays in synthesis-friendly way
            for (init_i = 0; init_i < NUM_BOXES; init_i = init_i + 1) begin
                scores[init_i] <= 0;
                for (init_j = 0; init_j < 5; init_j = init_j + 1) begin
                    boxes[init_i][init_j] <= 0;
                end
            end
            for (init_j = 0; init_j < 5; init_j = init_j + 1) begin
                temp_box[init_j] <= 0;
            end
        end else begin
            case (state)
                IDLE: begin
                    state <= LOAD_BOXES;
                    box_counter <= 0;
                    addr <= 0;
                    wait_counter <= 0;
                end
                
                LOAD_BOXES: begin
                    if (addr < TOTAL_ENTRIES) begin
                        addr <= addr + 1;
                    end else begin
                        // Wait for the last 2 data items to be processed
                        state <= WAIT_LAST_DATA;
                        wait_counter <= 0;
                    end
                    
                    // Store data using delayed_addr_2 (2-cycle delayed address)
                    // Only store when we have valid delayed address within range
                    if (delayed_addr_2 < TOTAL_ENTRIES && delayed_addr_2 != {ADDR_WIDTH{1'b1}}) begin
                        case (delayed_addr_2 % 5)
                            0: boxes[delayed_addr_2/5][0] <= data_out; // X
                            1: boxes[delayed_addr_2/5][1] <= data_out; // Y
                            2: boxes[delayed_addr_2/5][2] <= data_out; // W
                            3: boxes[delayed_addr_2/5][3] <= data_out; // H
                            4: begin
                                boxes[delayed_addr_2/5][4] <= data_out; // S
                                scores[delayed_addr_2/5] <= data_out;
                            end
                        endcase
                    end
                end
                
                WAIT_LAST_DATA: begin
                    // Wait 2 more cycles to get the last data from BRAM
                    if (wait_counter < 2) begin
                        wait_counter <= wait_counter + 1;
                        // Continue storing delayed data
                        if (delayed_addr_2 < TOTAL_ENTRIES && delayed_addr_2 != {ADDR_WIDTH{1'b1}}) begin
                            case (delayed_addr_2 % 5)
                                0: boxes[delayed_addr_2/5][0] <= data_out; // X
                                1: boxes[delayed_addr_2/5][1] <= data_out; // Y
                                2: boxes[delayed_addr_2/5][2] <= data_out; // W
                                3: boxes[delayed_addr_2/5][3] <= data_out; // H
                                4: begin
                                    boxes[delayed_addr_2/5][4] <= data_out; // S
                                    scores[delayed_addr_2/5] <= data_out;
                                end
                            endcase
                        end
                    end else begin
                        state <= SORT_INIT;
                        sort_passes <= 0;
                        sorting_done <= 0;
                    end
                end
                
                SORT_INIT: begin
                    box_counter <= 0;
                    state <= SORT_BOXES;
                end
                
                SORT_BOXES: begin
                    // Bubble sort implementation - one comparison per clock
                    if (box_counter < NUM_BOXES-1) begin
                        if (scores[box_counter] < scores[box_counter+1]) begin
                            // Store values in temp registers first
                            temp_score <= scores[box_counter];
                            temp_box[0] <= boxes[box_counter][0];
                            temp_box[1] <= boxes[box_counter][1];
                            temp_box[2] <= boxes[box_counter][2];
                            temp_box[3] <= boxes[box_counter][3];
                            temp_box[4] <= boxes[box_counter][4];
                            swap_needed <= 1'b1;
                            state <= SORT_SWAP;
                        end else begin
                            box_counter <= box_counter + 1;
                        end
                    end else begin
                        // One pass completed
                        sort_passes <= sort_passes + 1;
                        if (sort_passes < NUM_BOXES-1) begin
                            // Need more passes
                            state <= SORT_INIT;
                        end else begin
                            // Sorting complete
                            state <= INIT_NMS;
                            current_box <= 0;
                            suppressed_reg <= {NUM_BOXES{1'b0}};
                            valid <= {NUM_BOXES{1'b1}};
                        end
                    end
                end
                
                SORT_SWAP: begin
                    if (swap_needed) begin
                        // Perform the swap using temp registers
                        scores[box_counter] <= scores[box_counter+1];
                        scores[box_counter+1] <= temp_score;
                        boxes[box_counter][0] <= boxes[box_counter+1][0];
                        boxes[box_counter][1] <= boxes[box_counter+1][1];
                        boxes[box_counter][2] <= boxes[box_counter+1][2];
                        boxes[box_counter][3] <= boxes[box_counter+1][3];
                        boxes[box_counter][4] <= boxes[box_counter+1][4];
                        boxes[box_counter+1][0] <= temp_box[0];
                        boxes[box_counter+1][1] <= temp_box[1];
                        boxes[box_counter+1][2] <= temp_box[2];
                        boxes[box_counter+1][3] <= temp_box[3];
                        boxes[box_counter+1][4] <= temp_box[4];
                        swap_needed <= 1'b0;
                    end
                    box_counter <= box_counter + 1;
                    state <= SORT_BOXES;
                end
                
                INIT_NMS: begin
                    $display("DEBUG: INIT_NMS - current_box=%0d, NUM_BOXES=%0d, valid[%0d]=%b, suppressed[%0d]=%b", 
                             current_box, NUM_BOXES, current_box, 
                             (current_box < NUM_BOXES) ? valid[current_box] : 1'b0, 
                             current_box, 
                             (current_box < NUM_BOXES) ? suppressed_reg[current_box] : 1'b0);
                    if (current_box < NUM_BOXES) begin
                        if (valid[current_box] && !suppressed_reg[current_box]) begin
                            compare_box <= current_box + 1;
                            $display("DEBUG: Starting comparison for box %0d with boxes %0d onwards", current_box, current_box + 1);
                            state <= COMPARE_BOXES;
                        end else begin
                            $display("DEBUG: Skipping box %0d (valid=%b, suppressed=%b)", 
                                     current_box, valid[current_box], suppressed_reg[current_box]);
                            current_box <= current_box + 1;
                        end
                    end else begin
                        $display("DEBUG: NMS complete - all boxes processed");
                        state <= FINISH;
                    end
                end
                
                COMPARE_BOXES: begin
                    $display("DEBUG: COMPARE_BOXES - compare_box=%0d, NUM_BOXES=%0d, valid=%b, suppressed=%b", 
                             compare_box, NUM_BOXES, valid[compare_box], suppressed_reg[compare_box]);
                    if (compare_box < NUM_BOXES) begin
                        if (valid[compare_box] && !suppressed_reg[compare_box]) begin
                            $display("DEBUG: Calculating IOU between box %0d and box %0d", current_box, compare_box);
                            iou_start <= 1'b1;
                            state <= CALCULATE_IOU;
                        end else begin
                            $display("DEBUG: Skipping comparison with box %0d (valid=%b, suppressed=%b)", 
                                     compare_box, valid[compare_box], suppressed_reg[compare_box]);
                            compare_box <= compare_box + 1;
                            // Stay in COMPARE_BOXES state to check next box
                        end
                    end else begin
                        $display("DEBUG: Finished comparing box %0d with all others", current_box);
                        current_box <= current_box + 1;
                        state <= INIT_NMS;
                    end
                end
                
                CALCULATE_IOU: begin
                    iou_start <= 1'b0;
                    state <= UPDATE_VALID; // Combinational IOU - go directly to update
                end
                
                UPDATE_VALID: begin
                    $display("DEBUG: UPDATE_VALID - IOU=%0d, threshold=%0d, will suppress=%b", 
                             iou_result, IOU_THRESHOLD, (iou_result > IOU_THRESHOLD));
                    if (iou_result > IOU_THRESHOLD) begin
                        $display("DEBUG: Suppressing box %0d", compare_box);
                        suppressed_reg[compare_box] <= 1'b1;
                        valid[compare_box] <= 1'b0;
                    end else begin
                        $display("DEBUG: Keeping box %0d", compare_box);
                    end
                    compare_box <= compare_box + 1;
                    state <= COMPARE_BOXES;
                end
                
                FINISH: begin
                    // NMS complete
                    state <= FINISH;
                end
                
                default: begin
                    state <= IDLE;
                end
            endcase
        end
    end

    // Debug display for BRAM access
    always @(posedge clk) begin
        if (state == LOAD_BOXES || state == WAIT_LAST_DATA) begin
            if (delayed_addr_2 < TOTAL_ENTRIES && delayed_addr_2 != {ADDR_WIDTH{1'b1}}) begin
                $display("ADDR = %2d, Data = %6d", delayed_addr_2, data_out);
            end
        end
    end

    // Debug display for loaded boxes
    always @(posedge clk) begin
        if (state == SORT_INIT && sort_passes == 0) begin
            $display("=== Loaded Boxes (Before Sorting) ===");
            for (init_i = 0; init_i < NUM_BOXES; init_i = init_i + 1) begin
                $display("Box[%0d]: X=%0d, Y=%0d, W=%0d, H=%0d, S=%0d", 
                         init_i, boxes[init_i][0], boxes[init_i][1], boxes[init_i][2], boxes[init_i][3], boxes[init_i][4]);
            end
        end
    end

    // Debug display for sorted boxes
    always @(posedge clk) begin
        if (state == INIT_NMS && current_box == 0) begin
            $display("=== Sorted Boxes (By Score) ===");
            for (init_i = 0; init_i < NUM_BOXES; init_i = init_i + 1) begin
                $display("Box[%0d]: X=%0d, Y=%0d, W=%0d, H=%0d, S=%0d", 
                         init_i, boxes[init_i][0], boxes[init_i][1], boxes[init_i][2], boxes[init_i][3], boxes[init_i][4]);
            end
        end
    end

    // Debug display for NMS progress
    always @(posedge clk) begin
        if (state == INIT_NMS) begin
            $display("NMS: Processing current_box = %0d, valid = %b, suppressed = %b", 
                     current_box, valid, suppressed_reg);
        end
        if (state == COMPARE_BOXES) begin
            $display("NMS: Comparing current_box = %0d with compare_box = %0d", 
                     current_box, compare_box);
        end
    end

    // Debug display for IOU calculations
    always @(posedge clk) begin
        if (state == UPDATE_VALID) begin
            $display("IOU between Box[%0d] and Box[%0d] = %0d (threshold=%0d) -> %s", 
                     current_box, compare_box, iou_result, IOU_THRESHOLD,
                     (iou_result > IOU_THRESHOLD) ? "SUPPRESS" : "KEEP");
        end
    end

    // Debug display for final suppression results
    always @(posedge clk) begin
        if (state == FINISH) begin
            $display("=== NMS Results ===");
            $display("Suppressed mask: %b", suppressed_reg);
            for (init_i = 0; init_i < NUM_BOXES; init_i = init_i + 1) begin
                $display("Box[%0d]: %s (Score=%0d)", init_i, 
                         suppressed_reg[init_i] ? "SUPPRESSED" : "KEPT", scores[init_i]);
            end
            $display("NMS completed at time %0t", $time);
            $display("Suppression results: %b", suppressed_reg);
        end
    end

endmodule

```
</details>

## Testbench (Verilog)
<details>
  <summary>All necessary signals are provided in testbench</summary>

```c
`timescale 1ns / 1ps

module nms_top_tb();

    // Parameters
    parameter NUM_BOXES = 8;
    parameter CLOCK_PERIOD = 10; // 10ns = 100MHz

    // Inputs
    reg clk;
    reg reset;

    // Outputs
    wire [NUM_BOXES-1:0] suppressed;

    // Instantiate the Unit Under Test (UUT)
    nms_top #(
        .NUM_BOXES(NUM_BOXES),
        .IOU_THRESHOLD(18'd180)  // 0.7 in Q8.8 format
    ) uut (
        .clk(clk),
        .reset(reset),
        .suppressed(suppressed)
    );

    // Clock generation
    initial begin
        clk = 0;
        forever #(CLOCK_PERIOD/2) clk = ~clk;
    end

    // Reset and stimulus
    initial begin
        // Initialize and reset
        reset = 1;
        #100;
        reset = 0;
        
        $display("Starting NMS simulation...");
        $display("Design will automatically read from block RAM");
        

        
        $finish;
    end
endmodule
```
</details>

# Non Maximum Suppression using general IOU formula with pipelining

## IOU Calculation general Algorithm (Verilog) with pipelining

<details>
  <summary>This is the sub moudle for NMS UNIT</summary>

```c
module calculate_iou (
    input clk,
    input rst,
    input  [17:0] x1, y1, w1, h1,
    input  [17:0] x2, y2, w2, h2,
    output reg [17:0] iou,
    output reg valid
);

    // Stage 1: Calculate box boundaries and intersection coordinates
    reg [17:0] x1w1_s1, x2w2_s1, y1h1_s1, y2h2_s1;
    reg [17:0] x1_s1, y1_s1, x2_s1, y2_s1, w1_s1, h1_s1, w2_s1, h2_s1;
    always @(posedge clk) begin
        if (rst) begin
            x1w1_s1 <= 0; x2w2_s1 <= 0; y1h1_s1 <= 0; y2h2_s1 <= 0;
            x1_s1 <= 0; y1_s1 <= 0; x2_s1 <= 0; y2_s1 <= 0;
            w1_s1 <= 0; h1_s1 <= 0; w2_s1 <= 0; h2_s1 <= 0;
        end else begin
            x1w1_s1 <= x1 + w1;
            x2w2_s1 <= x2 + w2;
            y1h1_s1 <= y1 + h1;
            y2h2_s1 <= y2 + h2;
            x1_s1 <= x1; y1_s1 <= y1; x2_s1 <= x2; y2_s1 <= y2;
            w1_s1 <= w1; h1_s1 <= h1; w2_s1 <= w2; h2_s1 <= h2;
        end
    end

    // Stage 2: Calculate intersection rectangle and areas
    reg [17:0] inter_x1_s2, inter_y1_s2, inter_x2_s2, inter_y2_s2;
    reg [17:0] inter_w_s2, inter_h_s2;
    reg [35:0] area1_s2, area2_s2;
    always @(posedge clk) begin
        if (rst) begin
            inter_x1_s2 <= 0; inter_y1_s2 <= 0; inter_x2_s2 <= 0; inter_y2_s2 <= 0;
            inter_w_s2 <= 0; inter_h_s2 <= 0;
            area1_s2 <= 0; area2_s2 <= 0;
        end else begin
            inter_x1_s2 <= (x1_s1 > x2_s1) ? x1_s1 : x2_s1;
            inter_y1_s2 <= (y1_s1 > y2_s1) ? y1_s1 : y2_s1;
            inter_x2_s2 <= (x1w1_s1 < x2w2_s1) ? x1w1_s1 : x2w2_s1;
            inter_y2_s2 <= (y1h1_s1 < y2h2_s1) ? y1h1_s1 : y2h2_s1;
            inter_w_s2  <= ( ((x1w1_s1 < x2w2_s1) ? x1w1_s1 : x2w2_s1) > ((x1_s1 > x2_s1) ? x1_s1 : x2_s1) ) ?
                            (((x1w1_s1 < x2w2_s1) ? x1w1_s1 : x2w2_s1) - ((x1_s1 > x2_s1) ? x1_s1 : x2_s1)) : 0;
            inter_h_s2  <= ( ((y1h1_s1 < y2h2_s1) ? y1h1_s1 : y2h2_s1) > ((y1_s1 > y2_s1) ? y1_s1 : y2_s1) ) ?
                            (((y1h1_s1 < y2h2_s1) ? y1h1_s1 : y2h2_s1) - ((y1_s1 > y2_s1) ? y1_s1 : y2_s1)) : 0;
            area1_s2 <= w1_s1 * h1_s1;
            area2_s2 <= w2_s1 * h2_s1;
        end
    end

    // Stage 3: Calculate intersection area and union
    reg [35:0] area_inter_s3, area1_s3, area2_s3;
    reg [35:0] area_union_s3;
    always @(posedge clk) begin
        if (rst) begin
            area_inter_s3 <= 0; area1_s3 <= 0; area2_s3 <= 0; area_union_s3 <= 0;
        end else begin
            area_inter_s3 <= inter_w_s2 * inter_h_s2;
            area1_s3 <= area1_s2;
            area2_s3 <= area2_s2;
            area_union_s3 <= area1_s2 + area2_s2 - (inter_w_s2 * inter_h_s2);
        end
    end

    // Stage 4: Calculate scaled IoU and output
    reg [43:0] scaled_iou_s4;
    always @(posedge clk) begin
        if (rst) begin
            scaled_iou_s4 <= 0;
            iou <= 0;
            valid <= 0;
        end else begin
            scaled_iou_s4 <= (area_union_s3 != 0) ? (area_inter_s3 * 256) / area_union_s3 : 0;
            iou <= scaled_iou_s4[17:0];
            valid <= 1'b1; // Output is valid after 4 cycles
        end
    end

endmodule

```
</details>


## NMS UNIT (Verilog)
<details>
  <summary>Compelete top module for NMS</summary>

```c
module nms_top #(
    parameter NUM_BOXES = 11,                  // Number of boxes (configurable)
    parameter IOU_THRESHOLD = 18'd180,         // 0.7 in Q8.8 format
    parameter DATA_WIDTH = 18,
    parameter ADDR_WIDTH = $clog2(NUM_BOXES * 5), // Dynamic address width
    parameter TOTAL_ENTRIES = NUM_BOXES * 5    // Dynamic total entries
)(
    input wire clk,
    input wire reset,
    output wire [NUM_BOXES-1:0] suppressed
);

    // Block RAM signals
    wire [DATA_WIDTH-1:0] data_out;
    reg [ADDR_WIDTH-1:0] addr;
    reg [ADDR_WIDTH-1:0] delayed_addr_1, delayed_addr_2;
    
    // 2-cycle delay for BRAM output to track which address the data corresponds to
    always @(posedge clk) begin
        if (reset) begin
            delayed_addr_1 <= {ADDR_WIDTH{1'b1}}; // Invalid address initially
            delayed_addr_2 <= {ADDR_WIDTH{1'b1}}; // Invalid address initially
        end else begin
            delayed_addr_1 <= addr;
            delayed_addr_2 <= delayed_addr_1;
        end
    end

    // Block RAM instance
    blk_mem_gen_0 coe_rom (
        .clka(clk),
        .ena(1'b1),
        .wea(1'b0),
        .addra(addr),
        .dina({DATA_WIDTH{1'b0}}),
        .douta(data_out)
    );

    // Box storage registers - dynamically sized
    reg [DATA_WIDTH-1:0] boxes [0:NUM_BOXES-1][0:4]; // X, Y, W, H, S for each box
    reg [DATA_WIDTH-1:0] scores [0:NUM_BOXES-1];
    reg [NUM_BOXES-1:0] valid;
    reg [NUM_BOXES-1:0] suppressed_reg;
    
    assign suppressed = suppressed_reg;

    // FSM states
    localparam IDLE = 4'd0;
    localparam LOAD_BOXES = 4'd1;
    localparam WAIT_LAST_DATA = 4'd2;
    localparam SORT_INIT = 4'd3;
    localparam SORT_BOXES = 4'd4;
    localparam SORT_SWAP = 4'd5;
    localparam INIT_NMS = 4'd6;
    localparam COMPARE_BOXES = 4'd7;
    localparam CALCULATE_IOU = 4'd8;
    localparam UPDATE_VALID = 4'd9;
    localparam FINISH = 4'd10;
    
    reg [3:0] state;
    reg [ADDR_WIDTH-1:0] box_counter;
    reg [ADDR_WIDTH-1:0] current_box;
    reg [ADDR_WIDTH-1:0] compare_box;
    reg [2:0] wait_counter;
    reg [$clog2(NUM_BOXES):0] sort_passes; // Dynamic width for sort passes
    reg sorting_done;
    reg swap_needed;
    
    // IOU calculation signals
    wire [DATA_WIDTH-1:0] iou_result;
    reg iou_start;
    wire iou_done = 1'b1; // Combinational - always done
    
    // Temporary registers for swapping
    reg [DATA_WIDTH-1:0] temp_score;
    reg [DATA_WIDTH-1:0] temp_box [0:4];
    
    // IOU calculator (combinational)
    calculate_iou iou_calculator (
        .x1(boxes[current_box][0]), .y1(boxes[current_box][1]), 
        .w1(boxes[current_box][2]), .h1(boxes[current_box][3]),
        .x2(boxes[compare_box][0]), .y2(boxes[compare_box][1]), 
        .w2(boxes[compare_box][2]), .h2(boxes[compare_box][3]),
        .iou(iou_result)
    );

    // Synthesis-friendly initialization
    integer init_i, init_j;
    
    // Main FSM
    always @(posedge clk) begin
        if (reset) begin
            state <= IDLE;
            box_counter <= 0;
            current_box <= 0;
            compare_box <= 0;
            suppressed_reg <= {NUM_BOXES{1'b0}};
            valid <= {NUM_BOXES{1'b1}};
            addr <= 0;
            iou_start <= 0;
            wait_counter <= 0;
            sort_passes <= 0;
            sorting_done <= 0;
            swap_needed <= 0;
            temp_score <= 0;
            
            // Initialize arrays in synthesis-friendly way
            for (init_i = 0; init_i < NUM_BOXES; init_i = init_i + 1) begin
                scores[init_i] <= 0;
                for (init_j = 0; init_j < 5; init_j = init_j + 1) begin
                    boxes[init_i][init_j] <= 0;
                end
            end
            for (init_j = 0; init_j < 5; init_j = init_j + 1) begin
                temp_box[init_j] <= 0;
            end
        end else begin
            case (state)
                IDLE: begin
                    state <= LOAD_BOXES;
                    box_counter <= 0;
                    addr <= 0;
                    wait_counter <= 0;
                end
                
                LOAD_BOXES: begin
                    if (addr < TOTAL_ENTRIES) begin
                        addr <= addr + 1;
                    end else begin
                        // Wait for the last 2 data items to be processed
                        state <= WAIT_LAST_DATA;
                        wait_counter <= 0;
                    end
                    
                    // Store data using delayed_addr_2 (2-cycle delayed address)
                    // Only store when we have valid delayed address within range
                    if (delayed_addr_2 < TOTAL_ENTRIES && delayed_addr_2 != {ADDR_WIDTH{1'b1}}) begin
                        case (delayed_addr_2 % 5)
                            0: boxes[delayed_addr_2/5][0] <= data_out; // X
                            1: boxes[delayed_addr_2/5][1] <= data_out; // Y
                            2: boxes[delayed_addr_2/5][2] <= data_out; // W
                            3: boxes[delayed_addr_2/5][3] <= data_out; // H
                            4: begin
                                boxes[delayed_addr_2/5][4] <= data_out; // S
                                scores[delayed_addr_2/5] <= data_out;
                            end
                        endcase
                    end
                end
                
                WAIT_LAST_DATA: begin
                    // Wait 2 more cycles to get the last data from BRAM
                    if (wait_counter < 2) begin
                        wait_counter <= wait_counter + 1;
                        // Continue storing delayed data
                        if (delayed_addr_2 < TOTAL_ENTRIES && delayed_addr_2 != {ADDR_WIDTH{1'b1}}) begin
                            case (delayed_addr_2 % 5)
                                0: boxes[delayed_addr_2/5][0] <= data_out; // X
                                1: boxes[delayed_addr_2/5][1] <= data_out; // Y
                                2: boxes[delayed_addr_2/5][2] <= data_out; // W
                                3: boxes[delayed_addr_2/5][3] <= data_out; // H
                                4: begin
                                    boxes[delayed_addr_2/5][4] <= data_out; // S
                                    scores[delayed_addr_2/5] <= data_out;
                                end
                            endcase
                        end
                    end else begin
                        state <= SORT_INIT;
                        sort_passes <= 0;
                        sorting_done <= 0;
                    end
                end
                
                SORT_INIT: begin
                    box_counter <= 0;
                    state <= SORT_BOXES;
                end
                
                SORT_BOXES: begin
                    // Bubble sort implementation - one comparison per clock
                    if (box_counter < NUM_BOXES-1) begin
                        if (scores[box_counter] < scores[box_counter+1]) begin
                            // Store values in temp registers first
                            temp_score <= scores[box_counter];
                            temp_box[0] <= boxes[box_counter][0];
                            temp_box[1] <= boxes[box_counter][1];
                            temp_box[2] <= boxes[box_counter][2];
                            temp_box[3] <= boxes[box_counter][3];
                            temp_box[4] <= boxes[box_counter][4];
                            swap_needed <= 1'b1;
                            state <= SORT_SWAP;
                        end else begin
                            box_counter <= box_counter + 1;
                        end
                    end else begin
                        // One pass completed
                        sort_passes <= sort_passes + 1;
                        if (sort_passes < NUM_BOXES-1) begin
                            // Need more passes
                            state <= SORT_INIT;
                        end else begin
                            // Sorting complete
                            state <= INIT_NMS;
                            current_box <= 0;
                            suppressed_reg <= {NUM_BOXES{1'b0}};
                            valid <= {NUM_BOXES{1'b1}};
                        end
                    end
                end
                
                SORT_SWAP: begin
                    if (swap_needed) begin
                        // Perform the swap using temp registers
                        scores[box_counter] <= scores[box_counter+1];
                        scores[box_counter+1] <= temp_score;
                        boxes[box_counter][0] <= boxes[box_counter+1][0];
                        boxes[box_counter][1] <= boxes[box_counter+1][1];
                        boxes[box_counter][2] <= boxes[box_counter+1][2];
                        boxes[box_counter][3] <= boxes[box_counter+1][3];
                        boxes[box_counter][4] <= boxes[box_counter+1][4];
                        boxes[box_counter+1][0] <= temp_box[0];
                        boxes[box_counter+1][1] <= temp_box[1];
                        boxes[box_counter+1][2] <= temp_box[2];
                        boxes[box_counter+1][3] <= temp_box[3];
                        boxes[box_counter+1][4] <= temp_box[4];
                        swap_needed <= 1'b0;
                    end
                    box_counter <= box_counter + 1;
                    state <= SORT_BOXES;
                end
                
                INIT_NMS: begin
                    $display("DEBUG: INIT_NMS - current_box=%0d, NUM_BOXES=%0d, valid[%0d]=%b, suppressed[%0d]=%b", 
                             current_box, NUM_BOXES, current_box, 
                             (current_box < NUM_BOXES) ? valid[current_box] : 1'b0, 
                             current_box, 
                             (current_box < NUM_BOXES) ? suppressed_reg[current_box] : 1'b0);
                    if (current_box < NUM_BOXES) begin
                        if (valid[current_box] && !suppressed_reg[current_box]) begin
                            compare_box <= current_box + 1;
                            $display("DEBUG: Starting comparison for box %0d with boxes %0d onwards", current_box, current_box + 1);
                            state <= COMPARE_BOXES;
                        end else begin
                            $display("DEBUG: Skipping box %0d (valid=%b, suppressed=%b)", 
                                     current_box, valid[current_box], suppressed_reg[current_box]);
                            current_box <= current_box + 1;
                        end
                    end else begin
                        $display("DEBUG: NMS complete - all boxes processed");
                        state <= FINISH;
                    end
                end
                
                COMPARE_BOXES: begin
                    $display("DEBUG: COMPARE_BOXES - compare_box=%0d, NUM_BOXES=%0d, valid=%b, suppressed=%b", 
                             compare_box, NUM_BOXES, valid[compare_box], suppressed_reg[compare_box]);
                    if (compare_box < NUM_BOXES) begin
                        if (valid[compare_box] && !suppressed_reg[compare_box]) begin
                            $display("DEBUG: Calculating IOU between box %0d and box %0d", current_box, compare_box);
                            iou_start <= 1'b1;
                            state <= CALCULATE_IOU;
                        end else begin
                            $display("DEBUG: Skipping comparison with box %0d (valid=%b, suppressed=%b)", 
                                     compare_box, valid[compare_box], suppressed_reg[compare_box]);
                            compare_box <= compare_box + 1;
                            // Stay in COMPARE_BOXES state to check next box
                        end
                    end else begin
                        $display("DEBUG: Finished comparing box %0d with all others", current_box);
                        current_box <= current_box + 1;
                        state <= INIT_NMS;
                    end
                end
                
                CALCULATE_IOU: begin
                    iou_start <= 1'b0;
                    state <= UPDATE_VALID; // Combinational IOU - go directly to update
                end
                
                UPDATE_VALID: begin
                    $display("DEBUG: UPDATE_VALID - IOU=%0d, threshold=%0d, will suppress=%b", 
                             iou_result, IOU_THRESHOLD, (iou_result > IOU_THRESHOLD));
                    if (iou_result > IOU_THRESHOLD) begin
                        $display("DEBUG: Suppressing box %0d", compare_box);
                        suppressed_reg[compare_box] <= 1'b1;
                        valid[compare_box] <= 1'b0;
                    end else begin
                        $display("DEBUG: Keeping box %0d", compare_box);
                    end
                    compare_box <= compare_box + 1;
                    state <= COMPARE_BOXES;
                end
                
                FINISH: begin
                    // NMS complete
                    state <= FINISH;
                end
                
                default: begin
                    state <= IDLE;
                end
            endcase
        end
    end

    // Debug display for BRAM access
    always @(posedge clk) begin
        if (state == LOAD_BOXES || state == WAIT_LAST_DATA) begin
            if (delayed_addr_2 < TOTAL_ENTRIES && delayed_addr_2 != {ADDR_WIDTH{1'b1}}) begin
                $display("ADDR = %2d, Data = %6d", delayed_addr_2, data_out);
            end
        end
    end

    // Debug display for loaded boxes
    always @(posedge clk) begin
        if (state == SORT_INIT && sort_passes == 0) begin
            $display("=== Loaded Boxes (Before Sorting) ===");
            for (init_i = 0; init_i < NUM_BOXES; init_i = init_i + 1) begin
                $display("Box[%0d]: X=%0d, Y=%0d, W=%0d, H=%0d, S=%0d", 
                         init_i, boxes[init_i][0], boxes[init_i][1], boxes[init_i][2], boxes[init_i][3], boxes[init_i][4]);
            end
        end
    end

    // Debug display for sorted boxes
    always @(posedge clk) begin
        if (state == INIT_NMS && current_box == 0) begin
            $display("=== Sorted Boxes (By Score) ===");
            for (init_i = 0; init_i < NUM_BOXES; init_i = init_i + 1) begin
                $display("Box[%0d]: X=%0d, Y=%0d, W=%0d, H=%0d, S=%0d", 
                         init_i, boxes[init_i][0], boxes[init_i][1], boxes[init_i][2], boxes[init_i][3], boxes[init_i][4]);
            end
        end
    end

    // Debug display for NMS progress
    always @(posedge clk) begin
        if (state == INIT_NMS) begin
            $display("NMS: Processing current_box = %0d, valid = %b, suppressed = %b", 
                     current_box, valid, suppressed_reg);
        end
        if (state == COMPARE_BOXES) begin
            $display("NMS: Comparing current_box = %0d with compare_box = %0d", 
                     current_box, compare_box);
        end
    end

    // Debug display for IOU calculations
    always @(posedge clk) begin
        if (state == UPDATE_VALID) begin
            $display("IOU between Box[%0d] and Box[%0d] = %0d (threshold=%0d) -> %s", 
                     current_box, compare_box, iou_result, IOU_THRESHOLD,
                     (iou_result > IOU_THRESHOLD) ? "SUPPRESS" : "KEEP");
        end
    end

    // Debug display for final suppression results
    always @(posedge clk) begin
        if (state == FINISH) begin
            $display("=== NMS Results ===");
            $display("Suppressed mask: %b", suppressed_reg);
            for (init_i = 0; init_i < NUM_BOXES; init_i = init_i + 1) begin
                $display("Box[%0d]: %s (Score=%0d)", init_i, 
                         suppressed_reg[init_i] ? "SUPPRESSED" : "KEPT", scores[init_i]);
            end
            $display("NMS completed at time %0t", $time);
            $display("Suppression results: %b", suppressed_reg);
        end
    end

endmodule
```
</details>

## Testbench (Verilog)
<details>
  <summary>All necessary signals are provided in testbench</summary>

```c
`timescale 1ns / 1ps

module nms_top_tb();

    parameter NUM_BOXES = 8;
    parameter CLOCK_PERIOD = 10; // 10ns = 100MHz

    // Inputs
    reg clk;
    reg reset;

    // Outputs
    wire [NUM_BOXES-1:0] suppressed;

    // Instantiate the Unit Under Test (UUT)
    nms_top #(
        .NUM_BOXES(NUM_BOXES),
        .IOU_THRESHOLD(18'd180)  // 0.7 in Q8.8 format
    ) uut (
        .clk(clk),
        .reset(reset),
        .suppressed(suppressed)
    );

    // Clock generation
    initial begin
        clk = 0;
        forever #(CLOCK_PERIOD/2) clk = ~clk;
    end

    // Reset and stimulus
    initial begin
        // Initialize and reset
        reset = 1;
        #100;
        reset = 0;

        $display("Starting NMS simulation...");
        $display("Design will automatically read from block RAM (ensure BRAM is initialized in simulation)");

        // Wait for completion (state == FINISH)
        wait(uut.state == 4'd10); // FINISH state in your FSM
        #1000;

        // Display results
        $display("NMS completed at time %0t ns", $time);
        $display("Suppression results: %b", suppressed);

        $finish;
    end

    // Optionally, monitor the suppressed output for changes
    always @(suppressed) begin
        $display("Suppressed mask changed: %b at time %0t", suppressed, $time);
    end

endmodule
```
</details>


## calculate_iou (Method 1- verilog)
<details>
  <summary>Date:4 july 2025</summary>

```c
module calculate_iou (
    input  [17:0] x1, y1, w1, h1,
    input  [17:0] x2, y2, w2, h2,
    output reg [17:0] iou
);
    reg [17:0] x1w1, x2w2, y1h1, y2h2;
    reg [17:0] dx, dy;
    reg [17:0] inter, denom;
    reg [35:0] scaled_inter;
    
    always @(*) begin
        // Calculate rectangle boundaries
        x1w1 = x1 + w1;
        x2w2 = x2 + w2;
        y1h1 = y1 + h1;
        y2h2 = y2 + h2;
        
        // Calculate intersection dimensions
        dx = (x1w1 < x2w2 ? x1w1 : x2w2) - (x1 > x2 ? x1 : x2);
        dy = (y1h1 < y2h2 ? y1h1 : y2h2) - (y1 > y2 ? y1 : y2);
        
        // Calculate intersection and union (using your formula)
        inter = dx + dy;
        denom = (w1 + h1 > w2 + h2) ? (w1 + h1) : (w2 + h2);
        
        // Scale and calculate final IOU (keeping your formula)
        scaled_inter = (denom != 0) ? (inter << 8) : 0;
        iou = (denom != 0) ? scaled_inter / denom : 18'd0;
    end

endmodule
```


</details>



## calculate_iou (Method 2- verilog)
<details>
  <summary>Date:4 july 2025</summary>

```c
module calculate_iou (
    input wire clk,
    input wire reset,
    input wire start,
    input  [17:0] x1, y1, w1, h1,
    input  [17:0] x2, y2, w2, h2,
    output reg [17:0] iou,
    output reg valid
);

    // Pipeline registers - Stage 1
    reg [17:0] x1_r1, y1_r1, w1_r1, h1_r1;
    reg [17:0] x2_r1, y2_r1, w2_r1, h2_r1;
    reg start_r1;
    
    // Pipeline registers - Stage 2
    reg [17:0] x1w1_r2, x2w2_r2, y1h1_r2, y2h2_r2;
    reg [17:0] x1_r2, y1_r2, x2_r2, y2_r2;
    reg start_r2;
    
    // Pipeline registers - Stage 3
    reg [17:0] dx_r3, dy_r3, denom_r3;
    reg start_r3;
    
    // Pipeline registers - Stage 4
    reg [17:0] inter_r4;
    reg [17:0] denom_r4;
    reg start_r4;
    
    // Pipeline registers - Stage 5
    reg [35:0] scaled_inter_r5;
    reg [17:0] denom_r5;
    reg start_r5;
    
    // Intermediate wires for combinational logic
    wire [17:0] x1w1_next, x2w2_next, y1h1_next, y2h2_next;
    wire [17:0] max_x_next, min_x_next, max_y_next, min_y_next;
    wire [17:0] dx_next, dy_next;
    wire [17:0] inter_next, denom_next;
    wire [35:0] scaled_inter_next;
    wire [17:0] iou_next;
    
    // Stage 1: Input registration and boundary calculation
    assign x1w1_next = x1_r1 + w1_r1;
    assign x2w2_next = x2_r1 + w2_r1;
    assign y1h1_next = y1_r1 + h1_r1;
    assign y2h2_next = y2_r1 + h2_r1;
    
    // Stage 2: Min/Max calculation (broken into smaller comparisons)
    assign max_x_next = (x1_r2 > x2_r2) ? x1_r2 : x2_r2;
    assign min_x_next = (x1w1_r2 < x2w2_r2) ? x1w1_r2 : x2w2_r2;
    assign max_y_next = (y1_r2 > y2_r2) ? y1_r2 : y2_r2;
    assign min_y_next = (y1h1_r2 < y2h2_r2) ? y1h1_r2 : y2h2_r2;
    
    // Stage 3: Intersection dimension calculation
    assign dx_next = (min_x_next > max_x_next) ? (min_x_next - max_x_next) : 18'd0;
    assign dy_next = (min_y_next > max_y_next) ? (min_y_next - max_y_next) : 18'd0;
    assign denom_next = (w1_r1 + h1_r1 > w2_r1 + h2_r1) ? (w1_r1 + h1_r1) : (w2_r1 + h2_r1);
    
    // Stage 4: Intersection area calculation
    assign inter_next = dx_r3 + dy_r3;
    
    // Stage 5: Scaling
    assign scaled_inter_next = (denom_r4 != 0) ? (inter_r4 * 18'd256) : 36'd0;
    
    // Stage 6: Final division
    assign iou_next = (denom_r5 != 0) ? scaled_inter_r5[17:0] / denom_r5[7:0] : 18'd0;
    
    // Pipeline implementation
    always @(posedge clk) begin
        if (reset) begin
            // Reset all pipeline registers
            x1_r1 <= 18'd0; y1_r1 <= 18'd0; w1_r1 <= 18'd0; h1_r1 <= 18'd0;
            x2_r1 <= 18'd0; y2_r1 <= 18'd0; w2_r1 <= 18'd0; h2_r1 <= 18'd0;
            start_r1 <= 1'b0;
            
            x1w1_r2 <= 18'd0; x2w2_r2 <= 18'd0; y1h1_r2 <= 18'd0; y2h2_r2 <= 18'd0;
            x1_r2 <= 18'd0; y1_r2 <= 18'd0; x2_r2 <= 18'd0; y2_r2 <= 18'd0;
            start_r2 <= 1'b0;
            
            dx_r3 <= 18'd0; dy_r3 <= 18'd0; denom_r3 <= 18'd0;
            start_r3 <= 1'b0;
            
            inter_r4 <= 18'd0; denom_r4 <= 18'd0;
            start_r4 <= 1'b0;
            
            scaled_inter_r5 <= 36'd0; denom_r5 <= 18'd0;
            start_r5 <= 1'b0;
            
            iou <= 18'd0;
            valid <= 1'b0;
        end else begin
            // Stage 1: Input registration
            x1_r1 <= x1; y1_r1 <= y1; w1_r1 <= w1; h1_r1 <= h1;
            x2_r1 <= x2; y2_r1 <= y2; w2_r1 <= w2; h2_r1 <= h2;
            start_r1 <= start;
            
            // Stage 2: Boundary calculation
            x1w1_r2 <= x1w1_next;
            x2w2_r2 <= x2w2_next;
            y1h1_r2 <= y1h1_next;
            y2h2_r2 <= y2h2_next;
            x1_r2 <= x1_r1; y1_r2 <= y1_r1;
            x2_r2 <= x2_r1; y2_r2 <= y2_r1;
            start_r2 <= start_r1;
            
            // Stage 3: Min/Max and intersection dimensions
            dx_r3 <= dx_next;
            dy_r3 <= dy_next;
            denom_r3 <= denom_next;
            start_r3 <= start_r2;
            
            // Stage 4: Intersection area
            inter_r4 <= inter_next;
            denom_r4 <= denom_r3;
            start_r4 <= start_r3;
            
            // Stage 5: Scaling
            scaled_inter_r5 <= scaled_inter_next;
            denom_r5 <= denom_r4;
            start_r5 <= start_r4;
            
            // Stage 6: Final result
            iou <= iou_next;
            valid <= start_r5;
        end
    end

endmodule
```

</details>

## coe file (bram)
<details>
  <summary>Date:5 july 2025</summary>

```c
memory_initialization_radix=10;
memory_initialization_vector=
7245,83482,81152,42675,229,
45235,84966,63258,29926,136,
7296,83226,80512,43290,226,
7270,82970,80614,43443,213,
47206,84582,63462,28314,209,
6502,84275,69683,69248,72,
34534,84122,67610,41037,54,
47002,84275,62899,28774,175,
8346,83328,81306,43878,191,
6938,83558,80614,42445,214,
12544,84147,69018,63206,53;
```



