from datetime import datetime
from Train import diffusion, reverse_diffusion

def main():
    formatted_time = datetime.now().strftime("%m%d%H%M")
    train_type = 'case' # ctrl
    add_method = 'code'
    batch_size = 16
    num_time_steps = 1000
    num_epochs = 150000
    save_epoch = 50000
    generate_num = 1000
    add_epochs = 0
    checkpoint_path = None
    # GPU required: This code is optimized for GPU execution. Running it on CPU may result in significantly slower performance.
    device = 'cuda:0' # None
    diffusion(batch_size=batch_size, num_time_steps=num_time_steps, num_epochs=num_epochs, save_epoch=save_epoch,
              add_epochs=add_epochs,add_method=add_method, train_type=train_type, formatted_time=formatted_time,
              checkpoint_path=checkpoint_path, device = device)
    checkpoint_path = './/checkpoints/epoch_' + str(num_epochs) + '_step' + str(num_time_steps) + '_' + str(
        train_type) + '_' + str(add_method) + '_' + str(formatted_time) + '.pt'
    reverse_diffusion(num_time_steps = num_time_steps, train_type=train_type, num_epochs = num_epochs, generate_num = generate_num,
                     add_method = add_method, formatted_time = formatted_time, checkpoint_path = checkpoint_path)
    print("training finished")

if __name__ == '__main__':
    main()
