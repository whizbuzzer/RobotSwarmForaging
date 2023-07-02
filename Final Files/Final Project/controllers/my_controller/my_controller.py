"""my_controller controller."""

# You may need to import some classes of the controller module. Ex:
#  from controller import Robot, Motor, DistanceSensor
from controller import Robot, DistanceSensor, Motor, Emitter, Receiver, LED
import random as rd
#import numpy as np 

TIME_STEP = 64

# create the Robot instance.
robot = Robot()

# Initializing distance sensors:
ps = []
psNames = ['ps0', 'ps1', 'ps2', 'ps3', 'ps4', 'ps5','ps6', 'ps7']


for i in range(8):
    ps.append(robot.getDistanceSensor(psNames[i]))
    ps[i].enable(TIME_STEP)
        
emitter = robot.getEmitter('emitter')
emitter.setRange(0.17)
bot_emitter = robot.getEmitter('r2r_emitter')
receiver = robot.getReceiver('receiver')
bot_receiver = robot.getReceiver('r2r_receiver')
receiver.enable(TIME_STEP)
bot_receiver.enable(TIME_STEP)
bot_led_g = robot.getLED('led8')  # walker indicator
bot_led_y = robot.getLED('led9')  # beacon indicator

# receiver.setChannel(1)
#msg_flag = 0  # to prevent continuous printing
# print(receiver.getSamplingPeriod())

  
food = 0
nest_detected = 0
foodPheromone = 0
nestPheromone = 0
time = 0
    
# Initializing motors:
leftMotor = robot.getMotor('left wheel motor')
rightMotor = robot.getMotor('right wheel motor')
leftMotor.setPosition(float('inf'))
rightMotor.setPosition(float('inf'))
leftMotor.setVelocity(0.0)
rightMotor.setVelocity(0.0)

role = rd.random()  # for deciding walker or beacon
#print(role)

# Main loop:
# - perform simulation steps until Webots is stopping the controller
while robot.step(TIME_STEP) != -1:

    time = robot.getTime()
    bot_led_g.set(1)    
    
    
    # Read the sensors:
    psValues = []
    
    # position sensors:
    for i in range(8):
        psValues.append(ps[i].getValue())
        
    left_obstacle = psValues[5]>80.0 or psValues[6]>80.0 or psValues[7]>80.0
    right_obstacle = psValues[0]>80.0 or psValues[1]>80.0 or psValues[2]>80.0    
    
    # Using info on obstacle to actuate wheels:
    MAX_SPEED = 6.28  # ~2*pi rad/s
    
    # Initializing motor speeds at 50% of MAX_SPEED:
    leftSpeed = 0.5 * MAX_SPEED
    rightSpeed = 0.5 * MAX_SPEED
    
    # Robot will turn right:
    if left_obstacle:
        leftSpeed += 0.5 * MAX_SPEED
        rightSpeed -= 0.5 * MAX_SPEED
    
    # Robot will turn left:    
    elif right_obstacle:
        leftSpeed -= 0.5 * MAX_SPEED
        rightSpeed += 0.5 * MAX_SPEED
              
    # Writing actuator outputs:
    leftMotor.setVelocity(leftSpeed)
    rightMotor.setVelocity(rightSpeed)
        
    #pass
    
    #role = 0.4
    if time > 30: #== 40: # and time < 10.1:
        #print("inside this")
        #print(role)
    #for time in np.linspace(10,11):
        if role < 0.5:
            leftMotor.setVelocity(0*MAX_SPEED)
            rightMotor.setVelocity(0*MAX_SPEED)
    
    left_vel = leftMotor.getVelocity()
    right_vel = rightMotor.getVelocity()
    
    beacon = left_vel == 0 and right_vel == 0
    if beacon:
        bot_led_g.set(0)
        bot_led_y.set(1)
        beacon_indicate = str.encode("I am beacon")
        bot_emitter.send(beacon_indicate)
        #print(beacon_indicate)
    #if time > 100:
        #pass
        
        
        
    if receiver.getQueueLength() > 0:
        #print(receiver.getEmitterDirection())
        message = receiver.getData()
        message = message.decode("utf-8")
        # print(message)
            #print(receiver.getQueueLength())            
        if message == "food":
            food = 1
        #print(message)                     
        if message == "nest":
            nest_detected = 1
        receiver.nextPacket()
        
    if bot_receiver.getQueueLength() > 0:
        hola = bot_receiver.getData()
        hola = hola.decode("utf-8")
        #if hola == "I am beacon":
            #print(hola)        
              
    
    if food == 1 and nest_detected == 1:    
        place_food_message = str.encode("food received")
        emitter.send(place_food_message)
        #msg_flag = 0
        food = 0
    if food == 0:
        nest_detected = 0
    
    pass
# Enter here exit cleanup code.
