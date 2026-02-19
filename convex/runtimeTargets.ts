import { mutation, query } from "./_generated/server";
import { v } from "convex/values";

export const createJoinToken = mutation({
  args: {
    token: v.string(),
    tenant_id: v.string(),
    owner_user_id: v.string(),
    queue_name: v.string(),
    created_at: v.string(),
    expires_at: v.string(),
    node_name_hint: v.optional(v.string()),
    label: v.optional(v.string()),
  },
  handler: async (ctx, args) => {
    await ctx.db.insert("runtime_join_tokens", {
      ...args,
      used_at: undefined,
    });
    return { ok: true } as const;
  },
});

export const consumeJoinToken = mutation({
  args: {
    token: v.string(),
    consumed_at: v.string(),
  },
  handler: async (ctx, args) => {
    const tokenDoc = await ctx.db
      .query("runtime_join_tokens")
      .withIndex("by_token", (q) => q.eq("token", args.token))
      .first();
    if (!tokenDoc) {
      return { ok: false as const, reason: "not_found" as const };
    }
    if (tokenDoc.used_at) {
      return { ok: false as const, reason: "already_used" as const };
    }
    if (tokenDoc.expires_at <= args.consumed_at) {
      return { ok: false as const, reason: "expired" as const };
    }
    await ctx.db.patch(tokenDoc._id, { used_at: args.consumed_at });
    return {
      ok: true as const,
      join_token: {
        token: tokenDoc.token,
        tenant_id: tokenDoc.tenant_id,
        owner_user_id: tokenDoc.owner_user_id,
        queue_name: tokenDoc.queue_name,
        created_at: tokenDoc.created_at,
        expires_at: tokenDoc.expires_at,
        node_name_hint: tokenDoc.node_name_hint,
        label: tokenDoc.label,
      },
    };
  },
});

export const addRuntimeTarget = mutation({
  args: {
    tenant_id: v.string(),
    user_id: v.string(),
    target_id: v.string(),
    queue_name: v.string(),
    server_id: v.string(),
    registered_at: v.string(),
    name: v.optional(v.string()),
    metadata: v.optional(v.any()),
  },
  handler: async (ctx, args) => {
    await ctx.db.insert("runtime_targets", args);
    return { ok: true } as const;
  },
});

export const listTargets = query({
  args: {
    tenant_id: v.string(),
    user_id: v.string(),
  },
  handler: async (ctx, args) => {
    return ctx.db
      .query("runtime_targets")
      .withIndex("by_user", (q) =>
        q.eq("tenant_id", args.tenant_id).eq("user_id", args.user_id),
      )
      .collect();
  },
});

export const getTargetById = query({
  args: {
    target_id: v.string(),
  },
  handler: async (ctx, args) => {
    return ctx.db
      .query("runtime_targets")
      .withIndex("by_target_id", (q) => q.eq("target_id", args.target_id))
      .first();
  },
});

export const getActiveTarget = query({
  args: {
    tenant_id: v.string(),
    user_id: v.string(),
  },
  handler: async (ctx, args) => {
    const active = await ctx.db
      .query("runtime_active_targets")
      .withIndex("by_user", (q) =>
        q.eq("tenant_id", args.tenant_id).eq("user_id", args.user_id),
      )
      .first();
    if (!active) {
      return null;
    }
    return ctx.db
      .query("runtime_targets")
      .withIndex("by_target_id", (q) => q.eq("target_id", active.target_id))
      .first();
  },
});

export const setActiveTarget = mutation({
  args: {
    tenant_id: v.string(),
    user_id: v.string(),
    target_id: v.string(),
    updated_at: v.string(),
  },
  handler: async (ctx, args) => {
    const existing = await ctx.db
      .query("runtime_active_targets")
      .withIndex("by_user", (q) =>
        q.eq("tenant_id", args.tenant_id).eq("user_id", args.user_id),
      )
      .first();
    if (existing) {
      await ctx.db.patch(existing._id, {
        target_id: args.target_id,
        updated_at: args.updated_at,
      });
      return { updated: true } as const;
    }
    await ctx.db.insert("runtime_active_targets", args);
    return { updated: false } as const;
  },
});
