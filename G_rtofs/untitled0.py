#!/usr/bin/python
                
    
    

n = 3#int(input())
grid = ['--p','-m-','---'] 

bot_pos_i = [i for i, s in enumerate(grid) if 'm' in s][0]
bot_pos_j = grid[bot_pos_i].index('m')
prin_pos_i = [i for i, s in enumerate(grid) if 'p' in s][0]
prin_pos_j = grid[prin_pos_i].index('p')

bot_pos_i = n-1-bot_pos_i
prin_pos_i = n-1-prin_pos_i

while(1):
    if bot_pos_i<prin_pos_i:
        print("UP")
        bot_pos_i+=1
    elif bot_pos_i>prin_pos_i:
        print("DOWN")
        bot_pos_i-=1
    else:
         break

while(1):
    if bot_pos_j<prin_pos_j:
        print("RIGHT")
        bot_pos_j+=1
    elif bot_pos_j>prin_pos_j:
        print("LEFT")
        bot_pos_j-=1
    else:
         break
