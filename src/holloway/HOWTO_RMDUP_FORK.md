## RMDUP FORK

###CONCEPT
1. Decide what gets into the stack.
    1. tid/pos/mtid/mpos

##ALGORITHM
1. Fill the stack until the key is different.
2. For i in range(len(stack))
    1. For j in range(i + 1, len(stack))
        1. if(HD(bam_get_aux(stack[i], "BS"), bam_get_aux(stack[j])) < HD_THRESH)
            1. update_read(stack[j], stack[i]) // Update stack j with stack i to avoid iterating back through.
            2. bam_destroy1(stack[j])
            3. stack[j] = NULL
3. For i in range(len(stack))
    1. if(stack[i])
        1. bam_write1(stack[i])
        2. bam_destroy1(stack[i])
        3. stack[i] = NULL;
    2. stack->n = 0

