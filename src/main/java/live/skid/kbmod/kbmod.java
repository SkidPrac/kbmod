package live.skid.kbmod;

import net.jafama.*;
import org.bukkit.Bukkit;
import org.bukkit.ChatColor;
import org.bukkit.command.Command;
import org.bukkit.command.CommandExecutor;
import org.bukkit.command.CommandSender;
import org.bukkit.enchantments.Enchantment;
import org.bukkit.entity.Player;
import org.bukkit.event.EventHandler;
import org.bukkit.event.EventPriority;
import org.bukkit.event.Listener;
import org.bukkit.event.entity.EntityDamageByEntityEvent;
import org.bukkit.event.entity.EntityDamageEvent;
import org.bukkit.event.player.PlayerVelocityEvent;
import org.bukkit.plugin.java.JavaPlugin;
import org.bukkit.util.Vector;

import java.util.HashMap;

public class kbmod extends JavaPlugin implements Listener, CommandExecutor {
    double knockbackHorizontalFriction = 2.0D;
    double knockbackVerticalFriction = 2.0D;
    double knockbackHorizontal = 0.5D;
    double knockbackVertical = 0.4D;
    double knockbackVerticalLimit = 0.4D;
    double knockbackExtraHorizontal = 0.4D;
    double knockbackExtraVertical = 0.1D;
    HashMap<Player, Vector> playerKnockbackHashMap = new HashMap<>();

    @EventHandler(priority = EventPriority.LOWEST, ignoreCancelled = true)
    public void onPlayerVelocityEvent(PlayerVelocityEvent event) {
        if (!playerKnockbackHashMap.containsKey(event.getPlayer())) return;
        event.setVelocity(playerKnockbackHashMap.get(event.getPlayer()));
        playerKnockbackHashMap.remove(event.getPlayer());
    }

    @EventHandler(priority = EventPriority.MONITOR, ignoreCancelled = true)
    public void onEntityDamageEntity(EntityDamageByEntityEvent event) {
        if (event.getDamager() instanceof Player && event.getEntity() instanceof Player && !event.isCancelled() && event.getCause().equals(EntityDamageEvent.DamageCause.ENTITY_ATTACK)) {
            if (event.getDamage(EntityDamageEvent.DamageModifier.BLOCKING) != 0) {
                return;
            }

            if (!(event.getEntity() instanceof Player)) return;
            Player victim = (Player) event.getEntity();

            if (!(event.getDamager() instanceof Player)) return;
            if (!event.getCause().equals(EntityDamageEvent.DamageCause.ENTITY_ATTACK)) return;
            if (event.getDamage(EntityDamageEvent.DamageModifier.BLOCKING) != 0) return;

            Player attacker = (Player) event.getDamager();

            double d0 = attacker.getLocation().getX() - victim.getLocation().getX();
            double d1 = attacker.getLocation().getZ() - victim.getLocation().getZ();

            for (d1 = attacker.getLocation().getZ() - victim.getLocation().getZ();
                 d0 * d0 + d1 * d1 < 1.0E-4D; d1 = (FastMath.random() - FastMath.random()) * 0.01D)
                d0 = (FastMath.random() - FastMath.random()) * 0.01D;

            double magnitude = FastMath.sqrt(d0 * d0 + d1 * d1);

            Vector playerVelocity = victim.getVelocity();

            playerVelocity.setX((playerVelocity.getX() / 2) - (knockbackHorizontalFriction * knockbackHorizontal));
            playerVelocity.setY((playerVelocity.getY() / 2) - (knockbackVerticalFriction * knockbackVertical));
            playerVelocity.setZ((playerVelocity.getZ() / 2) - (knockbackHorizontalFriction * knockbackHorizontal));

            int i = attacker.getItemInHand().getEnchantmentLevel(Enchantment.KNOCKBACK);
            if (attacker.isSprinting()) ++i;

            if (playerVelocity.getY() > knockbackVerticalLimit)
                playerVelocity.setY(knockbackVerticalLimit);

            if (i > 0)
                playerVelocity.add(new Vector((-FastMath.sin(attacker.getLocation().getYaw() * 3.1415927F / 180.0F) *
                        (float) i * knockbackExtraHorizontal), knockbackExtraVertical,
                        FastMath.cos(attacker.getLocation().getYaw() * 3.1415927F / 180.0F) *
                                (float) i * knockbackExtraHorizontal));
            playerKnockbackHashMap.put(victim, playerVelocity);
        }
    }

    @Override
    public void onEnable() {
        Bukkit.getPluginManager().registerEvents(this, this);
        Bukkit.getScheduler().runTaskTimer(this, playerKnockbackHashMap::clear, 1, 1);
    }
}
