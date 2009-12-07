# 1 "i386.c"
 
 
 
 
 
 

 
 

 


 

 












 








        .globl _addll  
        .globl _subll  
        .globl _addllx  
        .globl _subllx  
        .globl _shiftl  
        .globl _shiftlr  
        .globl _bfffo  
        .globl _mulll  
        .globl _addmul  
        .globl _divll  
        .globl _mulmodll  
        .globl _overflow  
        .globl _hiremainder  

.text

	.align 2
_addll:  
        xorl    %edx,%edx           
        movl    4(%esp),%eax        
        addl    8(%esp),%eax        
        adcl    %edx,%edx           
        movl    %edx,_overflow      
	ret                         
	.align 2,0x90

	.align 2
_addllx:  
        xorl    %edx,%edx           
        xorl    %eax,%eax           
        subl    _overflow  ,%eax    
        movl    4(%esp),%eax        
        adcl    8(%esp),%eax        
        adcl    %edx,%edx           
        movl    %edx,_overflow      
	ret                         
	.align 2,0x90

	.align 2
_subll:  
        xorl    %edx,%edx           
        movl    4(%esp),%eax        
        subl    8(%esp),%eax        
        adcl    %edx,%edx           
        movl    %edx,_overflow      
	ret                         
	.align 2,0x90

	.align 2
_subllx:  
        xorl    %edx,%edx           
        xorl    %eax,%eax           
        subl    _overflow  ,%eax    
        movl    4(%esp),%eax        
        sbbl    8(%esp),%eax        
        adcl    %edx,%edx           
        movl    %edx,_overflow      
	ret                         
	.align 2,0x90

	.align 2
_shiftl:  
        movl    4(%esp),%eax         
        movb    8(%esp),%cl          
        xorl    %edx,%edx            
        shldl   %cl,  %eax,%edx       

        shll    %cl,%eax             
        movl    %edx,_hiremainder    
	ret                          
	.align 2,0x90

	.align 2
_shiftlr:  
        movl    4(%esp),%eax         
        movb    8(%esp),%cl          
        xorl    %edx,%edx            
        shrdl   %cl,  %eax,%edx       

        shrl    %cl,%eax             
        movl    %edx,_hiremainder    
	ret                          
	.align 2,0x90

# 136 "i386.c"

	.align 2
_bfffo:  
        movl    4(%esp),%eax         
        bsrl    %eax,%edx            
        movl    $31,%eax
        subl    %edx,%eax            
        ret                          
	.align 2,0x90


	.align 2
_mulll:  
        movl    4(%esp),%eax         
        mull    8(%esp)              
        movl    %edx,_hiremainder    
        ret                          
	.align 2,0x90

	.align 2
_addmul:  
        xorl    %ecx,%ecx            
        movl    4(%esp),%eax         
        mull    8(%esp)              
        addl    _hiremainder  ,%eax  
        adcl    %ecx,%edx
        movl    %edx,_hiremainder    
        ret                          
	.align 2,0x90

        .align 2
_divll:  
        movl    4(%esp),%eax         
        movl    _hiremainder  ,%edx  
        divl    8(%esp)              
        movl    %edx,_hiremainder    
        ret                          
	.align 2,0x90

	.align 2
_mulmodll:  
        movl    4(%esp),%eax         
        mull    8(%esp)              
        divl    12(%esp)             
        movl    %edx,%eax            
        ret
	.align 2,0x90

.comm _overflow  ,4
.comm _hiremainder  ,4
